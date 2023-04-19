function avasmsplot(BestSol,varargin)
%avasmsplot produces the augmented star plot and enables interactivity
%
%<a href="matlab: docsearchFS('avasmsplot')">Link to the help function</a>
%
% Required input arguments:
%
%      BestSol :  Best solutions. table.
%                 A table containing the details of the admissible
%                 solutions which have been found. We define a solution as
%                 admissible if the residuals pass the Durbin-Watson and
%                 Jarque-Bera tests, at the 10 per cent level.
%                 If no solution is found, than a 5 per cent threshold for both
%                 tests is used.
%                 Note that in input option critBestSol it is possible to
%                 set up different thresholds to define the admissible
%                 solutions.
%                 The rows of BestSol are ordered in a non increasing way
%                 using the p-value of the Durbin Watson test (rescaled by
%                 the value of $R^22$ and the number of units not declared as
%                 outliers).
%                 Colums 1-5 contain boolean information about the usage of
%                 options: orderR2, scail, trapezoid, rob,
%                 tyinitial. These 5 columns are ordered in non decreasing
%                 way depending on the frequency in which a particular
%                 option has been used in the set of admissiable solutions.
%                 For example, if option orderR2 is the only one
%                 which is always used in the set of admissible solutions,
%                 then column 1 is associated to orderR2.
%                 6th column contains the value of R2 (column name R2).
%                 7th column contains the p-value of Durbin Watson
%                 test (column name pvalDW).
%                 8th column contains the p-value of Jarque Bera test
%                 (column name pvalJB).
%                 9th column contains the number of units which have not
%                 been declared as outliers (column name nused).
%                 10th column is a cell which contains the residuals for
%                 the associated solution.
%                 11th column contains the struct out which is the
%                 output of the call to avas.
%                 12th column contains the numbers obtained by the product
%                 of p-value of Durbin Watson test, the values of $R^2$ and
%                 the number of units which have not been declared as
%                 outliers. The values of this column depend on the
%                 optional input argument solOrdering.
%      Data Types - table
%
% Optional input arguments:
%
%  maxBestSol :  maximum number of admissible
%               solutions to show in the augmented star plot.
%               Missing or positive integer.
%               if maxBestSol is missing all solutions inside input table
%               BestSol are shown in the augmented star plot.
%               If maxBestSol is (say) 3, just the first 3 solutions
%                are shown in in the augmented star plot.
%           Example - 'maxBestSol',5
%           Data Types - double
%
%   showBars  : show bars of labels. Boolean.  If showBars is true
%               the values of R2, fraction of units used, pvalue of DW test
%               and pval of normality test are shown with bars below each
%               star, else (default) these values are shows using a
%               textbox. More precisely, If $k$ is the total
%               number of units declared as outliers the four bars
%               represent respectively $R^2$, $(n − 2k)/n$, $p_{DW}$ and
%               $p_{JB}$. Note that we have $2k$ rather than $k$ because the
%               maximum number of outliers is $[n/2]$.
%           Example - 'showBars',true
%           Data Types - logical
%
%  addPolygons : polygons around the outside. Boolean.
%               if addPolygons is true (default) , polygons around the
%               outside of the radial lines are added to the augmented star
%               plot. If addPolygons is false just the radial lines from
%               the center are shows.
%           Example - 'addPolygons',false
%           Data Types - logical
%
%
%       tag     :    Personalized plot tag. String. String which identifies
%                   the handle of the plot which
%                   is about to be created. The default is to use tag
%                   'pl_augstarplot'. Note that if the program finds a plot which
%                   has a tag equal to the one specified by the user, then
%                   the output of the new plot overwrites the existing one
%                   in the same window else a new window is created.
%                   Example - 'tag','myplot'
%                   Data Types - char
%
%    databrush  :   interactive mouse brushing. Empty value, scalar or structure.
%                   If databrush is an empty value (default), no brushing
%                   is done. The activation of this option (databrush is a
%                   scalar or a struct) enables the user  to select a set
%                   of solutions in the current plot and to see them shown
%                   in the aceplot (i.e. a plot which enables to visualize
%                   the results of avas).
%                   DATABRUSH IS A SCALAR
%                   If databrush is a scalar the default selection tool is
%                   a circular brush and it is possible to brush only once
%                   (that is persist='').
%                   DATABRUSH IS A STRUCTURE.
%                   If databrush is a structure, it is possible to use all
%                   optional arguments of function selectdataFS.m and the
%                   following optional argument:
%                   - persist. Persist is an empty value or a scalar
%                     containing the strings 'on' or 'off'.
%                     The default value of persist is '', that is brushing
%                     is allowed only once.
%                     If persist is 'on' or 'off' brushing can be done as
%                     many time as the user requires.
%                     If persist='on' then the aceplot for the solution(s) currently brushed
%                     are added to those previously brushed.
%                     If persist='off' every time a new brush is performed
%                     the aceplot containing previous solutions is removed.
%                   Example - 'databrush',1
%                   Data Types - single | double | struct
%
% Output:
%
% See also: avasms.m, avas.m
%
% References:
%
%   Riani M. and Atkinson A.C. and Corbellini A. (2023), Robust Transformations
%   for Multiple Regression via Additivity and Variance Stabilization,
%   submitted.
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('avasmsplot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% Example from Wang and Murphy: avasmsplot with all default options.
    close all
    rng('default')
    seed=100;
    negstate=-30;
    n=200;
    X1 = mtR(n,0,seed)*2-1;
    X2 = mtR(n,0,negstate)*2-1;
    X3 = mtR(n,0,negstate)*2-1;
    X4 = mtR(n,0,negstate)*2-1;
    res=mtR(n,1,negstate);
    % Generate y
    y = log(4 + sin(3*X1) + abs(X2) + X3.^2 + X4 + .1*res );
    X = [X1 X2 X3 X4];
    y([121 80 34 188 137 110 79 86 1])=1.9+randn(9,1)*0.01;
    
    % Automatic model selection
    [VALtfin,Resarraychk]=avasms(y,X,'plots',0);
    % Show the best solutions
    avasmsplot(VALtfin);
%}

%{
    %% Example of option maxBestSol.
    close all
    % Example from Wang and Murphy: 
    rng('default')
    seed=100;
    negstate=-30;
    n=200;
    X1 = mtR(n,0,seed)*2-1;
    X2 = mtR(n,0,negstate)*2-1;
    X3 = mtR(n,0,negstate)*2-1;
    X4 = mtR(n,0,negstate)*2-1;
    res=mtR(n,1,negstate);
    % Generate y
    y = log(4 + sin(3*X1) + abs(X2) + X3.^2 + X4 + .1*res );
    X = [X1 X2 X3 X4];
    y([121 80 34 188 137 110 79 86 1])=1.9+randn(9,1)*0.01;
    
    % Automatic model selection
    [VALtfin,Resarraychk]=avasms(y,X,'plots',0);
    % Show just the 4 best solutions
    avasmsplot(VALtfin,'maxBestSol',4);
%}


%{
    % Example of option addPolygons.
    close all
    % Example from Wang and Murphy: 
    rng('default')
    seed=100;
    negstate=-30;
    n=200;
    X1 = mtR(n,0,seed)*2-1;
    X2 = mtR(n,0,negstate)*2-1;
    X3 = mtR(n,0,negstate)*2-1;
    X4 = mtR(n,0,negstate)*2-1;
    res=mtR(n,1,negstate);
    % Generate y
    y = log(4 + sin(3*X1) + abs(X2) + X3.^2 + X4 + .1*res );
    X = [X1 X2 X3 X4];
    y([121 80 34 188 137 110 79 86 1])=1.9+randn(9,1)*0.01;
    
    % Automatic model selection
    [VALtfin,Resarraychk]=avasms(y,X,'plots',0);
    % Show just the radial lines and not the polygons around the 
    % outside of the radial lines 
    addPolygons=false;
    avasmsplot(VALtfin,'maxBestSol',4,'addPolygons',false);
%}

%{
    % Example of option tag.
    close all
    % Example from Wang and Murphy: 
    rng('default')
    seed=100;
    negstate=-30;
    n=200;
    X1 = mtR(n,0,seed)*2-1;
    X2 = mtR(n,0,negstate)*2-1;
    X3 = mtR(n,0,negstate)*2-1;
    X4 = mtR(n,0,negstate)*2-1;
    res=mtR(n,1,negstate);
    % Generate y
    y = log(4 + sin(3*X1) + abs(X2) + X3.^2 + X4 + .1*res );
    X = [X1 X2 X3 X4];
    y([121 80 34 188 137 110 79 86 1])=1.9+randn(9,1)*0.01;
    
    % Automatic model selection
    [VALtfin,Resarraychk]=avasms(y,X,'plots',0);
    % Set a personalized tag
    avasmsplot(VALtfin,'tag','mytag');
%}

%{
    %Interactive_example,
    % databrush is passed as a scalar.
    % In this case it is possible to brush only once.
    close all
    % Example from Wang and Murphy: 
    rng('default')
    seed=100;
    negstate=-30;
    n=200;
    X1 = mtR(n,0,seed)*2-1;
    X2 = mtR(n,0,negstate)*2-1;
    X3 = mtR(n,0,negstate)*2-1;
    X4 = mtR(n,0,negstate)*2-1;
    res=mtR(n,1,negstate);
    % Generate y
    y = log(4 + sin(3*X1) + abs(X2) + X3.^2 + X4 + .1*res );
    X = [X1 X2 X3 X4];
    y([121 80 34 188 137 110 79 86 1])=1.9+randn(9,1)*0.01;
    
    % Automatic model selection
    [VALtfin,Resarraychk]=avasms(y,X,'plots',0);
    
     % databush is passed as a scalar
     % In this case it is possible to brush only once
    avasmsplot(VALtfin,'databrush',1)
%}

%{
    % Interactive_example,
    % databrush is passed as a struct.
    close all
    % Example from Wang and Murphy: 
    rng('default')
    seed=100;
    negstate=-30;
    n=200;
    X1 = mtR(n,0,seed)*2-1;
    X2 = mtR(n,0,negstate)*2-1;
    X3 = mtR(n,0,negstate)*2-1;
    X4 = mtR(n,0,negstate)*2-1;
    res=mtR(n,1,negstate);
    % Generate y
    y = log(4 + sin(3*X1) + abs(X2) + X3.^2 + X4 + .1*res );
    X = [X1 X2 X3 X4];
    y([121 80 34 188 137 110 79 86 1])=1.9+randn(9,1)*0.01;
    
    % Automatic model selection
    [VALtfin,Resarraychk]=avasms(y,X,'plots',0);
    
    % databrush is passed as a struct
    databrush=struct;
    databrush.persist='on';
    databrush.BrushShape='rect';
    databrush.BrushSize=0.2;
    avasmsplot(VALtfin,'databrush',databrush)
%}

%% Beginning of code

if nargin < 1
    error('FSDA:avasmsplot:missingInputs', ...
        'Missing input table containing required arguments to compute avasmsplot.')
end

maxBestSol=[];
databrush='';
tag='pl_augstarplot';
showBars=false;
addPolygons=true;
options=struct('maxBestSol',maxBestSol,'tag',tag,...
    'databrush',databrush,'showBars',showBars,'addPolygons',addPolygons);

if nargin > 1
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    % Check if number of supplied options is valid
    if length(varargin) ~= 2*length(UserOptions)
        error('FSDA:avasmsplot:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
    end
    % Check if user options are valid options
    chkoptions(options,UserOptions)

    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    maxBestSol=options.maxBestSol;
    tag=options.tag;
    databrush=options.databrush;
    showBars=options.showBars;
    addPolygons=options.addPolygons;
end

% the maximum number of solutions to show is equalù
% to the rows of BestSol
if isempty(maxBestSol)
    maxSol=size(BestSol,1);
else
    maxSol=min([size(BestSol,1),maxBestSol]);
end

rowlabs="R2="+string(num2str(BestSol{:,"R2"},3))...
    +" n="+string(num2str(BestSol{:,"nused"}))...
    + newline ...
    +" dw="+string(num2str(BestSol{:,"pvalDW"},2))...
    +" jb="+string(num2str(BestSol{:,"pvalJB"},2));
rowlabs=rowlabs(1:maxSol,:);

varlabs=BestSol.Properties.VariableNames(1:5);
VALtadj=BestSol{1:maxSol,1:5}.*BestSol{1:maxSol,"ord"};

if showBars==true
    testdata=BestSol(1:maxSol,:);
else
    testdata='';
end
% call augStarplot with options BestSols and addPolygons
centers=augStarplot(VALtadj(1:maxSol,:),rowlabs(1:maxSol,:),varlabs, ...
    'BestSol',testdata,'addPolygons',addPolygons);

set(gcf,'Tag',tag,'Name', 'Augmented star plot', 'NumberTitle', 'off')


%% Brush mode (call to function selectdataFS)
if ~isempty(databrush) || isstruct(databrush)
    if isstruct(databrush)

        fdatabrush=fieldnames(databrush);

        % persist option
        if isfield(databrush,'persist')
            persist=databrush.persist;
            databrush=rmfield(databrush,'persist');
            fdatabrush=fieldnames(databrush);
        else
            persist='';
        end

        % transform the input structure databrush into a cell array
        cv=[fdatabrush struct2cell(databrush)]';
        sele=[cv(:)' 'Ignore' {findobj(gcf,'tag','env')}];
        % add the FlagSize of the brushed points if it has not been
        % previously specified by the user
        d=find(strcmp('FlagSize',fdatabrush));
        if d==0
            sele=[sele 'FlagSize' '3'];
        end
    else
        sele={'selectionmode' 'Brush'};
        persist='';
    end

    % some local variables
    but=0; brushcum=[]; ij=1;

    plot1=gcf;

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
                disp('Unique solutions selected so far')
                disp(brushcum)
            end
        end

        % call to selectdataFS
        disp('Interactive part: select the center of one or more stars');
        % [pl,xcoo,ycoo]=selectdataFS('selectionmode','Brush','brushsize',brushWidth,'FlagColor','k');
        [pl,xcoo,ycoo]=selectdataFS(sele{:});

        % exit if the avasmsplot figure was closed before selection
        if isnumeric(pl) && ~isempty(pl) && (pl == -999)
            return
        end

        %% For each brushing operation, do the following:
        if ~isempty(cell2mat(pl))

            % From the stars which have been selected
            cand=unique(cell2mat([xcoo ycoo]),'rows');
            [~,soli]=intersect(centers,cand,'rows');
            outj=BestSol{soli,"Out"};
            solj=outj{:};
            aceplot(solj,'oneplot',true)

            solitxt=['Solution ' num2str(soli)];
            sgtitle([solitxt regexprep(rowlabs(soli,:),'\n','')])
            disp(solitxt)
            disp(rowlabs(soli,:))


            % brushcum = - the list of selected solution in all
            % iterations if
            %   persist=on
            % - the list of selected observations in the current iteration
            %   if persist=off
            if strcmp(persist,'on')
                brushcum=unique([brushcum; soli]);
            else
                brushcum=soli;
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
                % - the augmented star plot is highlighted again
                figure(plot1);
                % - a function to be executed on figure close is set
                set(gcf,'CloseRequestFcn',@closereqFS);

                % - and lay down the plots before continuing
                position(plot1);
                disp('Highlight the augmented star plot then: click on it to continue brushing or press a keyboard key to stop');
                ss=waitforbuttonpressFS;
                disp('------------------------');

                % After waitforbuttonpress: - the standard MATLAB function
                % to be executed on figure
                % close is recovered
                set(gcf,'CloseRequestFcn','closereq');
                Open_tyX = findobj(0, 'type', 'figure','tag','pl_tyX');
                Open_augstarplot = findobj(0, 'type', 'figure','tag',tag);

                if strcmp('off',persist)
                    if ~isempty(Open_tyX)
                        delete(Open_tyX);
                    end    % tyX plot is deleted
                end

                if isempty(Open_augstarplot)  % User closed the main brushing window
                    if ~isempty(Open_tyX)
                        delete(Open_tyX);
                    end    % tyX plot is deleted
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
            position(plot1);
        else
            disp('Wrong selection: Try again');
            disp('Select the center of a star');
            figure(plot1);
        end % for each brushing operation do ...
    end % close loop associated with but (loop brushing)
end % close options.databrush
end
%FScategory:VIS-Reg