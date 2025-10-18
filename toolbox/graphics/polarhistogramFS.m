function H = polarhistogramFS(Y, varargin)
%creates a polarhistogram
%
%<a href="matlab: docsearchFS('polarhistogramFS')">Link to the help function</a>
%
%  Required input arguments:
%
%           Y: Input data. 2D array or table.
%               n x p data matrix; n observations and p variables. Rows of
%               Y represent observations, and columns represent variables.
%               Rows of TBL represent observations, and columns
%               represent variables.
%                Data Types - single|double
%
%  Optional input arguments:
%
%   DataVars :  Table variables for which to compute polarhistograms.
%               Numeric vector or string array of cell array of character
%               vectors. For examples if 'DataVars' is [3 5], polar
%               histogram is done for variables 3 and 5. If 'DataVars' is
%               ["Name1" "Name4"] variable with these names inside the
%               input table are used
%               Example - 'DataVars',[2 4]
%               Data Types -  character vector | string array | cell array of character vectors | vector of positive integers | logical vector
%
%    groupvars : Grouping Variable.  Identifiers for the grouping variables in input array Y.
%               If groupvars is [] then the output refers to the overall
%               sample. It can be vector of length n or a string (char) or
%               a number identifying a particular column in the input table
%               Y.
%                   Example - 'groupvars',2
%                   Data Types - character vector | string array | cell array of character vectors | vector of positive integers | logical vector | []
% 
%                 
%     nbins    : number of bins or bin edges. Scalar or vector.
%                If nbins is a scalar, then we assume that it is referred
%                to the number of bins. Alternatively if nbins is a numeric
%                vector of length>1, we assume that  nbins(1) is the leading
%                edge of the first bin, and nbins(end) is the trailing edge
%                of the last bin elements. The elements of input vector y
%                are binned into nbins equally spaced containers if nbins
%                is a scalar or into length(nbins)-1 containers if nbins is
%                not a scalar.
%               Example - 'nbins',10
%               Data Types - numeric vector
%                Remark: note that it is possible to pass all the options
%                which are allowed inside polarhistogram. See the examples
%                below for further details.
%
%  Output:
%
%    H       :  array of graphic objects. Graphic object of size
%               length(groupvars)-by-length(DataVars).  In position i,j it
%               contains the handle for ith level of groupvars for variable
%               DataVars.
%
% See also polarhistogram
%
% References:
%
%   Tufte E.R. (1983), "The visual display of quantitative information",
%   Graphics Press, Cheshire.
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('polarhistogramFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Example with all default arguments.
    load citiesItaly
    polarhistogramFS(citiesItaly);
%}

%{
    %% Call with grouping variable.
    load citiesItaly
    zone=[repelem("N",46) repelem("CS",57)]';
    polarhistogramFS(citiesItaly,'groupvars',zone);
%}

%{
    % Call with grouping variable and DataVars.
    load citiesItaly
    zone=[repelem("N",46) repelem("CS",57)]';
    polarhistogramFS(citiesItaly,'DataVars',[2 5],'groupvars',zone);
%}

%{
    % Call with grouping variable, DataVars and nbins.
    load citiesItaly
    zone=[repelem("N",46) repelem("CS",57)]';
    polarhistogramFS(citiesItaly,'DataVars',1:4,'groupvars',zone,'nbins',10);
%}

%{
    % Use of FaceAlpha.
    load citiesItaly
    zone=[repelem("N",46) repelem("CS",57)]';
    polarhistogramFS(citiesItaly,'DataVars',1:4,'groupvars',zone,'nbins',10,'FaceAlpha',0.01);
%}

%{
    % Grouping variable with more than two levels.
    load citiesItaly
    zone=[repelem("N",46) repelem("C",21) repelem("S",36)]';
    polarhistogramFS(citiesItaly,'DataVars',1:4,'groupvars',zone,'nbins',10);
%}

%{
    % Example of  DisplayStyle 'stairs'.
    load citiesItaly
    zone=[repelem("N",46) repelem("C",21) repelem("S",36)]';
    polarhistogramFS(citiesItaly,'DataVars',1:4,'groupvars',zone,'nbins',10,'DisplayStyle','stairs');
%}

%{
    %% Example of  DisplayStyle 'stairs' with Linewidth.
    load citiesItaly
    zone=[repelem("N",46) repelem("C",21) repelem("S",36)]';
    polarhistogramFS(citiesItaly,'DataVars',1:4,'groupvars',zone,'nbins',10,'DisplayStyle','stairs','LineWidth',2);
%}

%{
    load citiesItaly
    % Example of nbins as a vector.
    % In this case nbins contains the edges
    polarhistogramFS(citiesItaly,'DataVars',1:2,'groupvars','','nbins',[1000:1000:5000 30000],'DisplayStyle','stairs','LineWidth',2);
%}

%{
    % Call with grouping variable inside the table.
    load citiesItaly
    zone=[repelem("N",46) repelem("CS",57)]';
    C=citiesItaly;
    C.zone=zone;
    polarhistogramFS(C,'DataVars',2:5,'groupvars',"zone");
%}


%% Beginning of code
DataVars=[];
nbins=[];
groupvars=[];

if nargin>1
    options=struct('DataVars',DataVars,'nbins',nbins,'groupvars',groupvars);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid.
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:polarhistogramFS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end

        % Check if user options are valid options.
        % aux.chkoptions(options,UserOptions)
    end

    % Write in structure 'options' the options chosen by the user.
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    DataVars=options.DataVars;
    nbins=options.nbins;
    groupvars=options.groupvars;


    % Check if plots is present inside varargin
    checkDataVars = strcmp(UserOptions,'DataVars')>0;
    if any(checkDataVars)==true
        fcheckOutputFormat=2*find(checkDataVars);
        varargin([fcheckOutputFormat-1 fcheckOutputFormat])=[];
        UserOptions=varargin(1:2:length(varargin));

    end

    % Check if plots is present inside varargin
    checkGroupars = strcmp(UserOptions,'groupvars')>0;
    if any(checkGroupars)==true
        fcheckOutputFormat=2*find(checkGroupars);
        varargin([fcheckOutputFormat-1 fcheckOutputFormat])=[];
        UserOptions=varargin(1:2:length(varargin));
    end

    % Check if nbins is present inside varargin
    checknbins = strcmp(UserOptions,'nbins')>0;
    if any(checknbins)==true
        fcheckOutputFormat=2*find(checknbins);
        varargin([fcheckOutputFormat-1 fcheckOutputFormat])=[];
    end

end

if isempty(nbins)
    nbins=6;
end

n=size(Y,1);


if isempty(groupvars)
    groupvars=string(groupvars);
elseif length(groupvars)==n % In this case grouping variable is external to the table
    groupvars=string(groupvars);
elseif  isnumeric(groupvars) | ischar(groupvars) % In this case grouping variable is inside the table or the array
    if ischar(groupvars)
        group=Y{:,groupvars};
        Y(:,groupvars)=[];
        groupvars=group;
    elseif isnumeric(groupvars) & ~isempty(groupvars)
        if istable(Y)
            group=Y{:,groupvars};
        else
            group=Y(:,groups);
        end
        Y(:,groupvars)=[];
        groupvars=group;
    else
        error('FSDA:polarhistogramFS:WrongInputOpt','Could not identify grouping variable')
    end
else
    error('FSDA:polarhistogramFS:WrongInputOpt','Could not identify grouping variable')
end

p=size(Y,2);
if isempty(DataVars)
    DataVars=1:p;
elseif iscellstr(DataVars) | isstring(DataVars) % DataVars is a string
    if istable(Y)
        DataVarsInput=string(DataVars);
        [~,DataVars]=ismember(DataVarsInput,Y.Properties.VariableNames);
        if min(DataVars)==0
            error('FSDA:polarhistogramFS:WrongInputOpt','Variable names not all found inside the input table')
        end
    else
        error('FSDA:polarhistogramFS:WrongInputOpt','Variable names specified as strings but input is not a table')
    end
else

end

if istable(Y)
    names=Y.Properties.VariableNames;
    Y=Y{:,:};
else
    names="V"+(1:p);
end

ugr=unique(groupvars);

tiledlayout("flow")

% Normalize in the interval (0 2*pi) (extremes excluded)
[Y02pi,C,S]=normalize(Y(:,DataVars),"range",[1e-10 2*pi-1e-10]);
names=names(DataVars);



H=gobjects(max(1,length(ugr)),length(DataVars));


for j=1:length(DataVars)
    nexttile;

    if isscalar(nbins)
        edgesPolar=linspace(0,2*pi,nbins+1);
    else
        minj=min(Y(:,DataVars(j))); maxj=max(Y(:,DataVars(j)));
        nbins1=nbins;
        nbins1(nbins<minj | nbins>maxj)=[]; 
        if numel(nbins) < 1
            error('FSDA:polarhistogramFS:WrongInputOpt', ...
                'Number of supplied options is invalid. Wrong edges, 2 or more edges are needed. None of the edges is inside the interval [min max]');
        end
        % Normalize edgesPolar inside [0 2*pi]
        edgesPolar=normalize([minj nbins1 maxj],"Range",[0 2*pi]);
        if length(edgesPolar) > 3
            edgesPolar=edgesPolar(2:(end-1));
        end
    end


    if all(groupvars=="")
        h=polarhistogram(Y02pi(:,j),edgesPolar,varargin{:});
        H(1,j)=h;
    else
        for i=1:length(ugr)
            h=polarhistogram(Y02pi(groupvars==ugr(i),j), edgesPolar,varargin{:});
            hold on
            H(i,j)=h;
        end
    end
    % edgesPolar,'FaceAlpha',FaceAlpha,'Normalization',Normalization,'DisplayStyle',DisplayStyle
    set(gca,"ThetaTick",rad2deg(edgesPolar));
    edgesOriginal=edgesPolar*S(:,j)+C(:,j); % /(2*pi)
    if edgesPolar(end)>=2*pi && edgesPolar(1)==0
        set(gca,'ThetaTickLabel',string(edgesOriginal(1:end-1)));
    else
        set(gca,'ThetaTickLabel',string(edgesOriginal));
    end
    title(names(j))
    if j==length(DataVars) &&  ~isempty(ugr)
        legend(ugr)
    end

end

end
%FScategory:VIS-Reg