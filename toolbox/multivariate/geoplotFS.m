function geoplotFS(Y,score,ShapeFile)
%geoplotFS calls geoplotAPP.mlapp to show an interactive geoplot with polygons
%
%<a href="matlab: docsearchFS('geoplotFS')">Link to the help function</a>
% 
%
% In the app it is also possible choose the variabile of Y to show, the
% type of plot (i.e. in latitude and longitude coordinates or in the
% original coordinates of ShapeFile), the direction of the colorbar.
%
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
%  ShapeFile  : name of ShapeFile or geospatial table containg shapes. 
%               Character or string or geotable.
%               Name of the ShapeFile containing the containing the
%               geometric details of the rows. The ShapeFile, which is
%               loaded using function readgeotable, must have n rows and
%               the n rows must have the same order of the n rows of Y.  The function which is used to show
%               the plot depends on the fact that the mapping toolbox is
%               installed or not. If the mapping toolbox is installed
%               function geoscatter is called and the plot uses latitudes
%               and longitudes. On the other hand, if the mapping toolbox
%               is not installed function mapshow is used and the plot is
%               shown in shape coordinates without projection.
%                   Example - 'ShapeFile','shapefileName'
%                    Data Types - char or string or geotable.
% 
%  Optional input arguments:
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
% Copyright 2008-2024.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('geoplotFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:


%{
    %% Example where second input argument is a numeric scalar.
    
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

%% Beginning of code
  geoplotAPP(Y,score,ShapeFile)
end
%FScategory:VIS-Mult