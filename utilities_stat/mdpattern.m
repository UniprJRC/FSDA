function [Mispat,tMisAndOut] = mdpattern(Y, varargin)
%mdpattern finds and plots missing data patterns
%
% This function is useful for investigating any structure of missing
% observations in the data. In specific case, the missing data pattern
% could be (nearly) monotone. Monotonicity can be used to simplify the
% imputation model. See Schafer (1997) for details. Also, the missing
% pattern could suggest which variables could potentially be useful for
% imputation of missing entries.
%
%<a href="matlab: docsearchFS('mdpattern')">Link to the help function</a>
%
%  Required input arguments:
%
%     Y : data matrix (2D array)  table or timetable containing $n$
%         observations on $v$ quantitative variables.
%               Data Types - matrix, table or timetable.
%
%  Optional input arguments:
%
%
%       Lc   :  Vector of column labels. Cell array of charaters or String array.
%               Lc contains the labels of the columns of the input
%               array Y. This option is unnecessary if Y is a
%               table, because in this case Lc=X.Properties.VariableNames;
%               Example - 'Lc',{'Y1' Y2' 'Y3' 'Y4'}
%               Data Types - cell array of characters or String.
%
%
%       plots : Plot on the screen. Boolean
%               If plots = true (default), a plot which displays missing
%               data patterns is dispalyed on the screen.
%               Top axis contains the name of the variables
%               Big circle means missing value; smaller filled dot
%               represents non missing value Left axis shows the number of
%               observations for each pattern. For example number 40 shows
%               that the associated patterns is repeated 40 times.
%               The sum of the numbers on the left axis is n the total
%               number of rows
%               Right axis counts the variables with missing values and it
%               is equal to the number of big circles in the
%               corresponding row.
%                 Example - 'plots',false
%                 Data Types - Boolean
%
%  dispresults :  Display results on the screen. Boolean.
%                 If dispresults is true it is possible to see on the
%                 screen the two output tables Xpat,tMisAndOut.
%                 The default value of dispresults is false.
%                 Example - 'dispresults',true
%                 Data Types - Boolean
%
%
%
% Output:
%
%         Mispat:  missing values pattern. table.
%                table with size (k+1)x(v+2), where k is the total number
%                of missing values patterns which are present in the data
%                matrix. The first k rows contain the patterns.
%                The last row contains n and then the total number of
%                missing values in each column. The first column contains
%                information about the number of observations for each pattern.
%                The columns of Mispat are sorted in non decreasing number of
%                outliers. The last column contains the number of variables
%                with missing values for each pattern.
%
%   tMisAndOut:  missing values and univariate outliers for each variable. .
%                The rows of this table are associated with the variables.
%                The columns are referred to a series of statistics.
%                More precisely:
%                Columns 1:4 contain mean and median, std deviation
%                and rescaled MAD (median absolute deviation).
%                Fifth column (Count_miss) contains the number of missing
%                values for each variable.
%                Sixth column (Percmiss) conatins the percentage of missing data
%                for each variable.
%                Seventh and eight column contain the number of outliers
%                respectively in the left and right tail of the
%                distribution. The criterion to decide whether a unit is
%                outlier is based on the boxplot concept, that is the
%                outliers are the units which are above x0.75+1.5*IQR or
%                below x0.25-1.5*IQR, where IQR is the interquartile range.
%
%
% See also balloonplot, bubblechart
%
% References:
%
% Schafer, J.L. (1997). "Analysis of Incomplete Multivariate Data". London: Chapman & Hall.
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('mdpattern')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    %% mdpattern with table input.
    % Load the nhanes data
    % The nhanes data is a dataset with 25 observations on the following 4 variables.
    % age, Age group (1=20-39, 2=40-59, 3=60+).
    % bmi, Body mass index (kg/m**2).
    % hyp, Hypertensive (1=no,2=yes).
    % chl, Total serum cholesterol (mg/dL).
    % namvar array of strings containing the names of the columns of X.
     namvar=["age"  "bmi" "hyp" "chl"];
     X=[1   NaN  NaN  NaN
     2 22.7   1 187
     1   NaN   1 187
     3   NaN  NaN  NaN
     1 20.4   1 113
     3   NaN  NaN 184
     1 22.5   1 118
     1 30.1   1 187
     2 22.0   1 238
     2   NaN  NaN  NaN
     1   NaN  NaN  NaN
     2   NaN  NaN  NaN
     3 21.7   1 206
     2 28.7   2 204
     1 29.6   1  NaN
     1   NaN  NaN  NaN
     3 27.2   2 284
     2 26.3   2 199
     1 35.3   1 218
     3 25.5   2  NaN
     1   NaN  NaN  NaN
     1 33.2   1 229
     1 27.5   1 131
     3 24.9   1  NaN
     2 27.4   1 186];
     Xtable=array2table(X,VariableNames=namvar);
     [Mispat,tMisAndOut]=mdpattern(Xtable);
%}

%{
    %% Example of the use of options dispresults and plots.
    % Load the nhanes data
    % The nhanes data is a dataset with 25 observations on the following 4 variables.
    % age, Age group (1=20-39, 2=40-59, 3=60+).
    % bmi, Body mass index (kg/m**2).
    % hyp, Hypertensive (1=no,2=yes).
    % chl, Total serum cholesterol (mg/dL).
    % namvar array of strings containing the names of the columns of X.
     namvar=["age"  "bmi" "hyp" "chl"];
     X=[1   NaN  NaN  NaN
     2 22.7   1 187
     1   NaN   1 187
     3   NaN  NaN  NaN
     1 20.4   1 113
     3   NaN  NaN 184
     1 22.5   1 118
     1 30.1   1 187
     2 22.0   1 238
     2   NaN  NaN  NaN
     1   NaN  NaN  NaN
     2   NaN  NaN  NaN
     3 21.7   1 206
     2 28.7   2 204
     1 29.6   1  NaN
     1   NaN  NaN  NaN
     3 27.2   2 284
     2 26.3   2 199
     1 35.3   1 218
     3 25.5   2  NaN
     1   NaN  NaN  NaN
     1 33.2   1 229
     1 27.5   1 131
     3 24.9   1  NaN
     2 27.4   1 186];
     Xtable=array2table(X,VariableNames=namvar);
     % Plot is not shown.
     plots=false;
     % option dispresults shows a detailed explanation of
     % the content of two output matrices in the command window.
     dispresults=true;
     [Mispat,tMisAndOut]=mdpattern(Xtable,'plots',false,'dispresults',dispresults);
%}

%{
     % Example of mdpattern with timetable input.
     TT = readtimetable('outages.csv');
     [A,B]=mdpattern(TT(:,["Loss" "Customers" ]))
%}

%{
    %% An example with 2 simulated patterns of missing values.
    close all
    n=10000;
    p=10;
    X=randn(n,p);
    % Create first missing  data pattern
    n1=300; n2=3;
    rowsWithMis=randsample(n,n1);
    colsWithMis=randsample(p,n2);
    X(rowsWithMis,colsWithMis)=NaN;
    % Create second missing  data pattern
    n1=120; n2=5;
    rowsWithMis=randsample(n,n1);
    colsWithMis=randsample(p,n2);
    X(rowsWithMis,colsWithMis)=NaN;
    mdpattern(X);
%}

%% Beginning of code

[n,p]=size(Y);

plots=true;
dispresults=false;
Lc='';
if nargin>1
    
    options=struct('Lc',Lc,'plots',plots,...
        'dispresults',dispresults);
    
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        UserOptions=varargin(1:2:length(varargin));
        if ~isempty(UserOptions)
            % Check if number of supplied options is valid
            if length(varargin) ~= 2*length(UserOptions)
                error('FSDA:mdpattern:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
            end
            % Check if user options are valid options
            chkoptions(options,UserOptions)
        end
        
        % Write in structure 'options' the options chosen by the user
        if nargin > 2
            for i=1:2:length(varargin)
                options.(varargin{i})=varargin{i+1};
            end
        end
        Lc  = options.Lc;
        plots=options.plots;
        dispresults=options.dispresults;
    end
end

if ~(istable(Y) || istimetable(Y))
    if isempty(Lc)
        % Lc=cellstr(num2str((1:p)'));
        Lc=strrep(strcat("Y",num2str((1:p)')),' ','');
    else
        % Check that the length of Lc is equal to p
        if length(Lc)~=p
            error('FSDA:mdpattern:WrongInputOpt','Wrong length of column labels');
        end
    end
    Yd=Y;
else
    if istimetable(Y)
        Y=timetable2table(Y,'ConvertRowTimes',false);
    end
    Yd=table2array(Y);
    Lc=Y.Properties.VariableNames;
end


misY=ismissing(Y);
numbermisEachcol=sum(misY,1);
[numbermisEachcolOrd,ordercols]=sort(numbermisEachcol,'ascend');
Ysorcols=misY(:,ordercols);
YsorcolsVarNames=Lc(ordercols);

Yboo=~Ysorcols;

[Yuni,~,ib]=unique(Yboo,'rows');
tab=tabulateFS(ib);

Yuniz=flipud([tab(:,2) Yuni sum(~Yuni,2)]);

Yfin=[Yuniz;[n numbermisEachcolOrd sum(numbermisEachcolOrd)]];


totpa=size(Yfin,1);
rownamXfin="Pattern"+(1:totpa)';
rownamXfin(end)="totPatOrMis";
vanamXfin=[{'NrowsWithPattern'} YsorcolsVarNames(:)' {'NvarWithMis'}];

Mispat=array2table(Yfin,'RowNames',rownamXfin,'VariableNames',vanamXfin);


MisAndOut=zeros(p,8);
quant=[0.25 0.75];
for j=1:p
    Yj=Yd(:,j);
    mj=mean(Yj,'omitnan');
    medianj=median(Yj,'omitnan');
    sigmaj=std(Yj,'omitnan');
    madj= 1.4826*mad(Yj,1);
    countMissingj=sum(isnan(Yj));
    percMissingj=100*countMissingj/n;
    Yjnomiss=Yj(~isnan(Yj));
    quart=quantile(Yjnomiss,quant);
    DI=quart(2)-quart(1);
    ThreshInf=quart(1)-1.5*DI;
    ThreshSup=quart(2)+1.5*DI;
    outSUPj=sum(Yjnomiss>ThreshSup);
    outINFj=sum(Yjnomiss<ThreshInf);
    
    tot=[mj medianj sigmaj madj countMissingj percMissingj outSUPj outINFj];
    MisAndOut(j,:)=tot;
    
end

namesVariablesMisAndOut={'Mean' 'Median'   'Stdev' 'MAD'  'Count_miss' ...
    'Perc_miss'    'outInf'    'outSup'};
tMisAndOut=array2table(MisAndOut,"RowNames",Lc,"VariableNames",namesVariablesMisAndOut);

%% dispresults and plots
Mispat11=num2str(Mispat{1,1});

if dispresults==true
    disp('Table which shows missing values patterns')
    disp(Mispat)
    disp('0 means missing value and 1 represents non missing value')
    disp('First column contains the number of observations for each pattern')
    disp(['For example number ' Mispat11 ' shows that the associated pattern is repeated ' Mispat11 ' times'])
    disp('The sum of the numbers in the first column is n, that is the total number of rows')
    disp('The last column shows the number of variables with missing values for that particular pattern')
    disp('------------------------')
    disp('Missing value and outlier report')
    disp(tMisAndOut)
    disp('Columns outInf and outSup contain the number of units which are')
    disp('above x0.75+1.5*IQR or below x0.25-1.5*IQR, where IQR is the interquartile range')
end
if plots==true
    fs=14;
    Ysel=Yfin(1:end-1,2:end-1);
    balloonplot(~Ysel);
    map = [ 0 0 0.3 ; % Personalized color map of blues
        0 0 0.4 ;
        0 0 0.5 ;
        0 0 0.6 ;
        0 0 0.8 ;
        0 0 1.0];
    colormap(map);
    colorbar('off');
    bubblesize([2 30]);
    set(gcf,'Name','Missing data pattern figure');
    h=gca;
    h.YTickLabel=flip(string(Yfin(1:end-1,1)));
    h.XTickLabel=string(Yfin(end,2:end-1));
    h.FontSize=fs;
    xlabel("Number of missing values for each variable");
    ylabel("Number of rows with a particular pattern");
    ax1=gca;

    ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none','XTick',[],'YTick',[]);
    set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
    % set the same Limits and Ticks on ax2 as on ax1;
    if verLessThan('matlab','9.11')
        set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
    else
        set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'),'TickDir','none');
    end
    set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', get(ax1, 'YTick'),'FontSize',fs);
    OppositeYTickLabels = string(flip(Yfin(1:end-1,end)));
    
    % Set the x-tick and y-tick  labels for the second axes
    set(ax2 , 'TickLength' ,[0 0] , 'XMinorTick', 'off', ...
        'XTickLabel', YsorcolsVarNames,...
        'YTickLabel',OppositeYTickLabels);
    
    % yyaxis right
    ylabel(ax2,'Number of variables with missing values')
    
    linkaxes; % this is to link the limits of the top and bottom axses
    
    % Plot explanation
    disp('Detailed explanation of the "Missing data pattern figure"')
    disp('Top axis contains the names of the variables.')
    disp('Big circle means missing value; smaller filled dot represents non missing value.')
    disp('Left axis shows the number of observations for each pattern')
    disp(['For example number ' Mispat11 ' shows that the associated pattern is repeated ' Mispat11 ' times.'])
    disp('The sum of the numbers on the left axis is n, the total number of rows.')
    disp('Right axis counts the variables with missing values and')
    disp('it is equal to the number of big circles in the corresponding row.')
    disp('The number of missing values for each variable is shown on the bottom axis.')
    
end

end
%FScategory:UTISTAT
