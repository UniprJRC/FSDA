function geoplotFS(Y,score,ShapeFile, varargin)
%geoplotFS calls geoplotAPP.mlapp to show an interactive geoplot with polygons
%
%<a href="matlab: docsearchFS('geoplotFS')">Link to the help function</a>
%
%
% In the app it is also possible choose the variable of Y to show, the
% type of plot (i.e. in latitude and longitude coordinates or in the
% original coordinates of ShapeFile), the direction of the colorbar.
% Note that this function requires that the "mapping toolbox" is installed.
%
%
%  Required input arguments:
%
% Y :           Input data. 2D array or table.
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
%      score  : variable to show or nxp matrix containing variables
%               to show.
%               Numeric scalar, of string or character identifying the
%               initial variable to use in the geoplot. Alternatively,
%               score can contain a nxp array for example with the p
%               principal components
%                   Example - 'ShapeFile','shapefileName'
%                    Data Types - char or string numeric scalar or nxp array.
%
%  ShapeFile  : name of ShapeFile or geotable containing shapes.
%               Character or string or geotable.
%               Name of the ShapeFile containing the containing the
%               geometric details of the rows. The ShapeFile, which is
%               loaded using function readgeotable, must have n rows and
%               the n rows must have the same order of the n rows of Y.
%               Remark: note that this option can be used just if the
%               "Mapping toolbox" is installed.
%                   Example - 'ShapeFile','shapefileName'
%                    Data Types - char or string or geotable.
%
%  Optional input arguments:
%
%  bsb     :    list of units forming the initial subset. Vector. If bsb is
%               not specified all unis are considered
%                   Example - 'bsb',1:20
%                    Data Types - single|double
%
% Output:
%
%
%
% See also: biplotFS
%
% References:
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('geoplotFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    %% Example where second input argument is a numeric scalar.
    a=struct2table(ver);
    MappingInstalled=any(string(a{:,1})=="Mapping Toolbox");
    if MappingInstalled ==true
        load citiesItaly2024.mat
        X=citiesItaly2024;
        ShapeFile=X.Properties.UserData{1};
        % Show the geoplot of the third variable of table X
        geoplotFS(X,3,ShapeFile)
    else
       disp('This function requires the "mapping toolbox" is installed') 
    end
%}

%{
    %% Example where second input argument is a string.
    load citiesItaly2024.mat
    X=citiesItaly2024;
    ShapeFile=X.Properties.UserData{1};
    geoplotFS(X,"QualLif",ShapeFile)
%}


%{
    %% Example where second input argument is a character.
    load citiesItaly2024.mat
    X=citiesItaly2024;
    ShapeFile=X.Properties.UserData{1};
    geoplotFS(X,'SpendingA',ShapeFile)
%}


%{
    %% Example where second input argument is a matrix.
    load citiesItaly2024.mat
    X=citiesItaly2024;
    out=pcaFS(X,'biplot',0,'dispresults',0,'plots',0);
    ShapeFile=X.Properties.UserData{1};
    % Color based on PCs
    geoplotFS(X,out.score,ShapeFile)
%}

%{
    %% Example of use of option bsb.
    load citiesItaly2024.mat
    X=citiesItaly2024;
    out=pcaFS(X,'biplot',0,'dispresults',0,'plots',0);
    ShapeFile=X.Properties.UserData{1};
    % Color based on PCs and use of option bsb
    % Select the provinces of Emilia-Romagna
    bsb=[33:40 99];
    geoplotFS(X,out.score,ShapeFile,'bsb',bsb)
%}


%% Beginning of code
bsb=1:size(Y,1);

if nargin>3
    options=struct('bsb',bsb);


    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)


        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:geoplotFS:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end

        % Check if all the specified optional arguments were present
        % in structure options
        % Remark: the nocheck option has already been dealt by routine
        % chkinputR
        inpchk=isfield(options,UserOptions);
        WrongOptions=UserOptions(inpchk==0);
        if ~isempty(WrongOptions)
            disp(strcat('Non existent user option found->', char(WrongOptions{:})))
            error('FSDA:geoplotFS:NonExistInputOpt','In total %d non-existent user options found.', length(WrongOptions));
        end
    end


    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end

    bsb=options.bsb;
    
end

geoplotAPP(Y,score,ShapeFile,bsb)

end
%FScategory:VIS-Mult