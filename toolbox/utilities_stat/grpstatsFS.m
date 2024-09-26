function [statTable]=grpstatsFS(X, groupvars, whichstats, varargin)
%grpstatsFS calls grpstats and reshapes the output in a much better way
%
%<a href="matlab: docsearchFS('grpstatsFS')">Link to the help function</a>
%
%   grpstatsFS calls grpstats but shows the output in much better way. The
%   output of grpstatsFS is a table with a number of rows equal to the
%   number of variables for which statistics are computed. The number of
%   columns of this table is equal to the number of statistics which are
%   computed. In presence of a grouping variable the number of rows of the
%   output table remains the same but the number of columns is equal to the
%   number of statistics times the number of groups.
%
%  Required input arguments:
%
%            X: Input data. Table. Table containing n observations on v variables.
%               Rows of Y represent observations, and columns
%               represent variables.
%
%    groupvars: grouping variable.
%               Identifiers for the grouping variables in input X.
%               If groupvars is [] than the output refers to the overall
%               sample. For additional information on groupvars see the
%               help of grpstats.
%                   Example - 'groupvars',[ones(50,1); 2*ones(50,1)]
%                   Data Types - character vector | string array | cell array of character vectors | vector of positive integers | logical vector | []
%
%  whichstats: Types of summary statistics.
%              Name of the statistics which have to be computed
%              For additional information on whichstats see the help of
%              grpstats. If whichstats is empty of is not specifies the
%              summary statistics which are computed are ["mean" "median"
%              "std" "MAD" "skewness" "kurtosis"]
%               Example - 'groupvars',"group"
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
%  Output:
%
% statTable : table with p rows. The rows are referred to the variables of the input table
%             and the columns to the requested statistics.
%             The number of columns of statTable is equal to
%             the number of required statistcs times
%             the number of groups.
%
%
% See also: grpstats
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
    %% grpstats with just one input argument.
    % Load a table
    load citiesItaly.mat
    % Compute mean, median,  std, MAD, skewness, kurtosis
    % for the 7 variables of the input table citiesItaly
    TBL=grpstatsFS(citiesItaly);
    disp(TBL)
%}

%{
    %% grpstats with second input the grouping variable.
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
    %% Example of call to grpstats with personalized statistics.
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
    %%  Example of call to grpstatsFS to create conf int for the mean
    load citiesItaly.mat
    % Note that in this case meanci has in output two columns
    stats={"meanci" 'mean'};
    % Confidence interval and the sample means
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
    %  Example of the use of option DataVars.
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

%% Beginning of code
if nargin<2
    groupvars=[];
end

% Normalize MAD is one of the statistics which is computed
mads=@(x)median(abs(x-median(x)))/norminv(0.75);

if nargin<3 || isempty(whichstats)
    whichstats={"@mean" "@median" "@std" mads "@skewness" "@kurtosis"};
    nomiStat=["mean" "median" "std" "MAD" "skewness" "kurtosis"];
else

    if iscell(whichstats)
        nomiStat=cellfun(@char,whichstats,'UniformOutput',false);
    else
        nomiStat=whichstats;
    end
    nomiStat=strrep(nomiStat,"@","");

    fnd=find(nomiStat=="meanci");
    if ~isempty(fnd)
        nomiStat=[nomiStat(1:(fnd-1)) "meanCIinf" "meanCIsup" nomiStat(fnd+1:end)];
    end
end

if ~isempty(varargin)
    UserOptions=varargin{1:2:length(varargin)};

    % Check if DataVars is present inside varargin
    checkDataVars = strcmp(UserOptions,'DataVars')>0;
    if checkDataVars==true
        lmsval = varargin{2*find(checkDataVars)};
        vnames=X.Properties.VariableNames(lmsval);
    else
        vnames=X.Properties.VariableNames;
    end
else
    if istable(X)
        vnames=X.Properties.VariableNames;
    else
        error('FSDA:grpstatsFS:WrongInp','grpstatsFS just supports a table in input')
        % vnames="X"+string(1:size(X,2));
    end
end
if ~isempty(groupvars)
    vnames=setdiff(vnames,groupvars,'stable');
end
p=length(vnames);
tabTutti=grpstats(X,groupvars,whichstats,varargin{:});
ngroups=size(tabTutti,1);
lstats=length(nomiStat);
if ngroups>1
    nomivar=nomiStat'+tabTutti{:,1}';
    nomivar=nomivar(:);

    statArray=zeros(p,lstats*ngroups);

    for j=1:ngroups
        statArray(:,(lstats*(j-1)+1):(lstats*j))=reshape(tabTutti{j,3:end},lstats,p)';
    end
    statTable=array2table(statArray,"RowNames",vnames,"VariableNames",nomivar);
else
    statArray=reshape(tabTutti{1,2:end},lstats,p)';
    statTable=array2table(statArray,"RowNames",vnames,"VariableNames",nomiStat);
end
end
