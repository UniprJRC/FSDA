function [statTable]=grpstatsFS(TBL, groupvars, whichstats, varargin)
%grpstatsFS calls grpstats and reshapes the output in a much better way
%
%<a href="matlab: docsearchFS('grpstatsFS')">Link to the help function</a>
%
%   grpstatsFS calls grpstats, but shows the output in much better way. The
%   output of grpstatsFS is a table with a number of rows equal to the
%   number of variables for which statistics are computed. The number of
%   columns of this table is equal to the number of statistics which are
%   computed. By default, the statistics which are computed are the two
%   (non robust and robust) indexes of location, (mean and median) the two
%   (non robust and robust) indexes of spread (standard deviation and
%   scaled MAD) and the two (non robust and robust) indexes of skewness.
%   The robust index of skewness is the medcouple. The scaled MAD is
%   defined as 1.4826(med|x-med(x)|).
%   In presence of a grouping variable the number of rows of the
%   output table remains the same, but the number of columns is equal to the
%   number of statistics multiplied by the number of groups. With option
%   'OutputFormat' 'nested' which is similar to option  'OutputFormat' in
%   function pivot, it is possible to display the results using nested
%   tables.
%
%  Required input arguments:
%
%           TBL: Input data. Table. Table containing n observations on p variables.
%               Rows of TBL represent observations, and columns
%               represent variables. If it necessary to compute the
%               statistics for subgroups, TBL must include at least one
%               grouping variable, which you specify using groupvars.
%
%    groupvars: grouping variable.
%               Identifiers for the grouping variables in input TBL.
%               If groupvars is [] then the output refers to the overall
%               sample. For additional information on groupvars see the
%               help of grpstats. For example
%                   Example - 'groupvars',2
%                   Data Types - character vector | string array | cell array of character vectors | vector of positive integers | logical vector | []
%
%  whichstats: Types of summary statistics.
%              Name of the statistics which have to be computed.
%              For additional information on whichstats see the help of
%              grpstats. If whichstats is empty or it is not specified, the
%              summary statistics which are computed are ["mean" "median"
%              "std" "MAD" "skewness" "medcouple"];
%               Example - ["mean" "std"]
%               Data Types -  character vector | string array | function handle | cell array of character vectors or function handles.
%
%  Optional input arguments:
%
%      Alpha : Significance level. Scalar in [0 1).
%              Significance level for confidence and prediction intervals.
%              For additional information on Alpha see the help of
%              grpstats.
%               Example - 'Alpha',0.01
%               Data Types -  double
%
%   DataVars :  Table variables for which to compute summary statistics.
%               For additional information on Alpha see the help of
%               grpstats.
%               Example - 'DataVars',[2 4]
%               Data Types -  character vector | string array | cell array of character vectors | vector of positive integers | logical vector
%
%    VarNames : Variable names for output table. cell array or characters
%               or string array.
%               Note that the length of VarNames must be equal to the
%               number of statistics which are computed. Variable (column)
%               names for the output table statTable, specified as a string
%               array or a cell array of character vectors. By default,
%               grpstatsFS removes the @ if it is present in the name of
%               the statistic. In presence of a grouping variable
%               grpstatsFS appends the name corresponding to each
%               category of the groups.
%               Example - 'VarNames',["location" "robust location"];
%               Data Types -  character vector | string array | cell array of character vectors | vector of positive integers | logical vector
%
% OutputFormat : Column hierarchy output format. 'flat' (default) | 'nested'
%                In presence of one or more grouping variable, this option
%                specifies whether the label of each statistic must be
%                concatenated with the categories of the grouping variables
%                or the overall output table must contain nested tables
%                each referred to a different statistic.
%               Example - 'OutputFormat','nested';
%               Data Types -  character | string
%
%  Output:
%
% statTable : table with p rows.
%             Table containing summary statistics for the table input TBL.
%             The rows are referred to the variables of the input table
%             and the columns to the requested statistics.
%             The number of columns of statTable is equal to
%             the number of requested statistics multiplied by
%             the number of groups.
%
%
% See also: grpstats, medcouple
%
% References:
%
%
%
% Copyright 2008-2024.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('grpstatsFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

%{
    %% grpstatsFS with just one input argument.
    % Load a table
    load citiesItaly.mat
    % Compute mean, median,  std, MAD, skewness and medcouple
    % for the 7 variables of the input table citiesItaly
    TBL=grpstatsFS(citiesItaly);
    disp(TBL)
%}

%{
    %% grpstatsFS with second input the grouping variable.
    load citiesItaly.mat
    % The first 46 rows are referred to provinces located in northern Italy and
    % the remaining in centre-south Italy.
    zone=[repelem("N",46) repelem("CS",57)]';
    % Add zone to citiesItaly
    citiesItaly.zone=zone;
    TBL=grpstatsFS(citiesItaly,"zone");
    % First 6 columns are referred to province of the north.
    % The remaning columns to the other provinces.
    disp(TBL)
%}

%{
    %% Example of call to grpstatsFS with personalized statistics.
    % The second element is empty, that is there is no grouping variable.
    load citiesItaly.mat
    % Just compute mean and median
    stats=["mean" "median"];
    TBL=grpstatsFS(citiesItaly,[],stats);
    disp(TBL)
%}

%{
    %%  Example of call to grpstatsFS with function handles.
    load citiesItaly.mat
    % multiple types of summary statistics specified as a cell
    stats={"mean","std",@skewness};
    TBL=grpstatsFS(citiesItaly,[],stats);
    disp(TBL)
%}

%{
    %%  Example of call to grpstatsFS to create conf int for the mean.
    load citiesItaly.mat
    % Note that in this case meanci has in output two columns
    stats={"meanci" 'mean'};
    % Confidence interval for the sample means
    TBL=grpstatsFS(citiesItaly,[],stats);
    disp(TBL)
%}

%{
    %%  Example of call to grpstatsFS to create conf int for the mean with groups.
    load citiesItaly.mat
    % Note that in this case meanci has in output two columns
    stats={"meanci" 'mean'};
    % The first 46 rows are referred to provinces located in northern Italy and
    % the remaining in centre-south Italy.
    zone=[repelem("N",46) repelem("CS",57)]';
    % Add zone to citiesItaly
    citiesItaly.zone=zone;
    % Confidence interval for the sample means separated for the 2 groups
    TBL=grpstatsFS(citiesItaly,"zone",stats);
    disp(TBL)
%}

%{
    %%  Example of use of option Alpha.
    load citiesItaly.mat
    % Note that in this case meanci has in output two columns
    stats={"meanci" 'mean'};
    % The first 46 rows are referred to provinces located in northern Italy and
    % the remaining in centre-south Italy.
    zone=[repelem("N",46) repelem("CS",57)]';
    % Add zone to citiesItaly
    citiesItaly.zone=zone;
    % 99 per cent confidence intervals for the sample means separated for the 2 groups
    TBL=grpstatsFS(citiesItaly,"zone",stats,'Alpha',0.01);
    disp(TBL)
%}

%{
    %%  Example of the use of option DataVars.
    load citiesItaly.mat
    % Note that in this case meanci has in output two columns
    stats={"meanci" 'mean'};
    % The first 46 rows are referred to provinces located in northern Italy and
    % the remaining in centre-south Italy.
    zone=[repelem("N",46) repelem("CS",57)]';
    % Add zone to citiesItaly
    citiesItaly.zone=zone;
    % Confidence interval for the sample means separated for the 2 groups
    TBL=grpstatsFS(citiesItaly,"zone",stats,'DataVars',["addedval" "unemploy"]);
    disp(TBL)
%}

%{
    %%  Example of the use of option DataVars with VarNames.
    load citiesItaly.mat
    % Note that in this case meanci has in output two columns
    stats={"median" 'mean'};
    TBL=grpstatsFS(citiesItaly,[],stats, ...
        'DataVars',[1 2 5],'VarNames', ...
        ["Robust location" "Non robust location"]);
    disp(TBL)
%}

%{
    %%  Example of the use of option DataVars with VarNames and grouping variable.
    load citiesItaly.mat
    % Note that in this case meanci has in output two columns
    stats={"meanci" 'mean'};
    % The first 46 rows are referred to provinces located in northern Italy and
    % the remaining in centre-south Italy.
    zone=[repelem("N",46) repelem("CS",57)]';
    % Add zone to citiesItaly
    citiesItaly.zone=zone;
    % Confidence interval for the sample means separated for the 2 groups
    TBL=grpstatsFS(citiesItaly,"zone",stats, ...
        'DataVars',["addedval" "unemploy"],'VarNames', ...
        ["Mean: lower confidence interval" "Mean: upper confidence interval" "Sample mean"]);
    disp(TBL)
%}

%{
    %%  Example of the use of option OutputFormat.
    load citiesItaly.mat
    % Two measures of location and dispersion
    stats={'mean' 'median' 'std' @(x)1.4826*mad(x,1)};
    % The first 46 rows are referred to provinces located in northern Italy and
    % the remaining in centre-south Italy.
    zone=[repelem("N",46) repelem("CS",57)]';
    % Add zone to citiesItaly
    citiesItaly.zone=zone;
    % Requested statistics for the 2 groups using nested tables
    TBL=grpstatsFS(citiesItaly,"zone",stats,'OutputFormat','nested');
    format bank
    disp(TBL)
%}

%% Beginning of code
if nargin<2
    groupvars=[];
else
    groupvars=string(groupvars);
end


% Normalize MAD is one of the statistics which is computed
mads=@(x)median(abs(x-median(x)))/norminv(0.75);

% Option which controls whether to show the results using or not nested
% tables
OutputFormat='flat';

if nargin<3 || isempty(whichstats)
    whichstats={"@mean" "@median" "@std" mads "@skewness" "@medcouple"};
    nomiStat=["mean" "median" "std" "MAD" "skewness" "medcouple"];
    predefined=true;
else
    predefined=false;
end

if nargin>=3

    [varargin{:}] = convertStringsToChars(varargin{:});

    if predefined ==false
        if iscell(whichstats)
            nomiStat=cellfun(@char,whichstats,'UniformOutput',false);
        else
            nomiStat=whichstats;
        end
        nomiStat=strrep(nomiStat,"@","");
    end

    fnd=find(nomiStat=="meanci");
    if ~isempty(fnd)
        nomiStat=[nomiStat(1:(fnd-1)) "meanCIinf" "meanCIsup" nomiStat(fnd+1:end)];
    end

end

if ~isempty(varargin)
    UserOptions=varargin(1:2:length(varargin));

    % Check if DataVars is present inside varargin
    checkDataVars = strcmp(UserOptions,'DataVars')>0;
    if any(checkDataVars==true)
        lmsval = varargin{2*find(checkDataVars)};
        vnames=TBL.Properties.VariableNames(lmsval);
    else
        vnames=TBL.Properties.VariableNames;
    end

    % Check if VarNames is present inside varargin
    checkVarNames = strcmp(UserOptions,'VarNames')>0;
    if any(checkVarNames)==true
        fvarNames=2*find(checkVarNames);
        lmsval = varargin{fvarNames};
        nomiStat=lmsval;
        varargin([fvarNames-1 fvarNames])=[];
    end

    % Check if OutputFormat is present inside varargin
    checkOutputFormat = strcmp(UserOptions,'OutputFormat')>0;
    if any(checkOutputFormat)==true
        fcheckOutputFormat=2*find(checkOutputFormat);
        lmsval = varargin{fcheckOutputFormat};
        OutputFormat=lmsval;
        varargin([fcheckOutputFormat-1 fcheckOutputFormat])=[];
    end


else
    if istable(TBL) || istimetable(TBL)
        vnames=TBL.Properties.VariableNames;
    else
        error('FSDA:grpstatsFS:WrongInp','grpstatsFS just supports a table or a timetable in input')
    end
end
if ~isempty(groupvars)
    vnames=setdiff(vnames,groupvars,'stable');
end
p=length(vnames);
tabTutti=grpstats(TBL,groupvars,whichstats,varargin{:});
lgroupvars=length(groupvars);
ngroups=size(tabTutti,1);

if ischar(nomiStat)
    nomiStat=string(nomiStat);
end

lstats=length(nomiStat);
if ngroups>1

    str=string(tabTutti{:,1:lgroupvars});

    str= join(str  ,"_",2);
    nomivar=nomiStat'+str'; %  tabTutti{:,1:lgroupsvars}';
    nomivar=nomivar(:);

    % nomivarNested is used just if OutputFormat is nested
    nomivarNested=nomiStat'+strings(1,length(str));
    nomivarNested=nomivarNested(:);

    statArray=zeros(p,lstats*ngroups);

    for j=1:ngroups
        statArray(:,(lstats*(j-1)+1):(lstats*j))=reshape(tabTutti{j,(2+lgroupvars):end},lstats,p)';
    end

    if strcmp(OutputFormat,'flat')==true
        statTable=array2table(statArray,"RowNames",vnames,"VariableNames",nomivar);
    else
        CE=cell(lstats,1);
        for i=1:lstats
            boo=strcmp(nomivarNested,nomiStat(i));
            nomiboo=nomivar(boo);
            nomiboo1=strrep(nomiboo,nomiStat(i),"");
            CE{i}=array2table(statArray(:,boo),"VariableNames",nomiboo1);
        end
        statTable=table(CE{:},'RowNames',vnames,'VariableNames',nomiStat);
    end

else
    statArray=reshape(tabTutti{1,2:end},lstats,p)';
    statTable=array2table(statArray,"RowNames",vnames,"VariableNames",nomiStat);
end
end
%FScategory:UTISTAT