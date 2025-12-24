function ax = dumbbellPlot(X, varargin)
%dumbbellPlot creates a dumbbell chart comparing two sets of values
%
%
%<a href="matlab: docsearchFS('dumbbellPlot')">Link to the help function</a>
%
% A dumbbell plot is a data visualization used to compare two values per
% category. Each category is represented by two points connected by a line,
% resembling a dumbbell. The two points usually indicate before–after, A vs
% B, or two time periods. The distance between the points highlights the
% magnitude of change or difference. Colors or shapes can distinguish the
% two values clearly. It is especially effective when comparing many
% categories without clutter. Dumbbell plots emphasize direction of change
% better than bar charts. They work well with ordered categories to show
% trends. They require a common numeric scale for meaningful comparison.
% Overall, dumbbell plots provide a clean and intuitive comparison tool.
%
%  Required input arguments:
%
%            X: Input data. Vector or table with one or more columns.
%               Data matrix containing n observations
%                 Data Types -  2D array or table 
%
%  Optional input arguments:
%
%  orientation  : orientation of the plot. String or char.
%                Admissible values are 'horizontal' (default) or 'vertical'
%               Example - 'Orientation','vertical'
%               Data Types - string or char
%
%   labelX1    : Custom legend label for the first set of data. Char or string. 
%               The Default label is "X1"  if X1 is numeric vector or the corresponding
%               table name if X1 is a table.
%               Example - 'labelX1','mylabelBefore'
%               Data Types - char or string
%
%   labelX2    : Custom legend label for the second set of data. Char or string. 
%               The Default label is "X2"  if X2 is numeric vector or the corresponding
%               table name if X2 is a table.
%               Example - 'labelX2','mylabelAfter'
%               Data Types - char or string
%
%   plotType   : Type of plot layout. Char or String.
%               Determines whether to create a single plot or a side by side
%               plot. Admissible values are 'single' (default) or 'double'
%               Example - 'plotType', 'double'
%               Data Types - char or string
%
%   Title      : plot title(s). string or char or cell array.
%               Title(s) for the chart(s). The number of input titles, must 
%               be equal to the number of plots that are created. Default 
%               value for double plots is ['Year 1' 'Year 2']
%               Example - 'Title', '2025 products revenue'
%               Data Types - string or char or cell array
%
%   YLabels    : Custom tick labels. string or char or cell array.
%              Custom labels for each category (row) in the plot. Must have
%              the same length as the number of rows. If not provided,
%              default names will be used ('Row 1', 'Row 2',...)
%              Example - 'YLabels', ["Product A", "Product B", "Product C"]
%              Data Types - string | char | cell array
%
%   Color      : Color of the dots. string or char or 2x3 numeric array.
%              Specifies the colors for the two sets of dots. Can be either
%              a built-in palette ('default', 'colorblind', 'ruby_jade','cherry_sky', 'red_blue')
%               or a 2x3 array of RGP triplets where each row defines the
%               color of a set of dots.
%              Example - 'Color', [0.88, 0.30, 0.30; 0.20, 0.42, 0.85]
%              Data Types - string | char | array of 2 valid MATLAB colours
%             
%   MarkerSize : Size of the dots. numeric scalar or vector.
%              If scalar, all markers use the same size. If vector, it must
%              have the same length as the number of data pairs (rows), and
%              each pair of dots will use the corresponding size value.
%              Example - 'MarkerSize', 200
%              Data Types - numeric scalar or vector
%
%   LineWidth  : Width of the bars of the plot. numeric scalar.
%               Controls the thickness of the lines (bars) connecting each
%               pair of dots in the dumbbell plot.
%               Example - 'LineWidth', 2.5
%               Data Types - Double
%
%   TextSize   : Font size of data value labels. numeric scalar.
%              Controls the font size of the numeric labels displayed
%              near or inside the dots.
%              Example - 'TextSize', 14
%              Data Types - double
%
%   TextInside : Position of data value labels. Logical scalar.
%              When true, numeric labels are displayed inside the marker dots.
%              if false, labels are positioned outside the markers. 
%              Example - 'TextInside', true
%              Data Types - logical
%
% AxesFontSize : Font size of axis labels and ticks. Numeric scalar.
%              Controls the font size of the axis labels.
%              Example - 'AxesFontSize', 15
%              Data Types - double
%
%   ColorDist  : Color mapping based on difference between points. String or char.
%              Controls whether and how the connecting lines are colored based
%              on the difference between the two values. 
%              Admissible values are:'false' (default) 'directional','magnitude','robust'
%              Example - 'ColorDist', 'directional'
%              Data Types - string | char
%
% Background : Style of background of the plot. String or char.
%              Controls whether alternating background bands are displayed
%              behind the plot. Admissible values: 'none' (default), 'bands'
%              Example - 'Background', 'bands'
%              Data Types - string | char
%
%
%
%  Output:
%
%      ax :    handles to the charts. Graphic handle or Vector of graphic handles.
%              If the plot is single, a single ax handle is return; else, ax
%              is a vector containing the graphics handles for all the
%              plots
%
%
% See also funnelchart
%
% References:
%
%   Cleveland, W. S. (1984). Graphical Methods for Data Presentation, 
%   "The American Statistician", Vol. 38, pp. 270–280
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%
%<a href="matlab: docsearchFS('dumbbellPlot')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:
%{
    %% Example data for a 10-row dumbbell chart.
    categories = "Item " + (1:10)';                 
    beforeVals  = [52; 47; 63; 58; 44; 71; 39; 66; 55; 49];
    afterVals   = [61; 50; 68; 62; 48; 78; 45; 70; 60; 54]; 
    T = table(beforeVals, afterVals, ...
        'VariableNames', {'Before','After'},'RowNames',categories);
    dumbbellPlot(T);
%}

%{
    %% Example double vertical plot.
    % Patient health metrics: systolic BP comparison across two time periods
    patients = {'Patient A', 'Patient B', 'Patient C', 'Patient D'};
    
    % Period 1: January measurements
    jan_morning = [120, 135, 128, 142];
    jan_evening = [128, 145, 135, 150];
    
    % Period 2: June measurements (after lifestyle changes)
    jun_morning = [115, 125, 120, 135];
    jun_evening = [122, 132, 128, 142];
    figure;
    dumbbellPlot(jan_morning, jan_evening, jun_morning, jun_evening, ...
        'plotType', 'double','orientation', 'vertical','labelX1', 'Morning', ...
      'labelX2', 'Evening','Title', {'January BP', 'June BP'},'YLabels', patients);
%}

%{
    %% Example of MarkerSize being used to show weights.
    products = {'Product A', 'Product B', 'Product C', 'Product D', 'Product E'};
    before_campaign = [65, 45, 78, 52, 88];
    after_campaign = [82, 68, 95, 71, 92];

    % Marker sizes based on revenue importance (arbitrary weights)
    importance = [300, 200, 150, 90, 450];

    figure;
    dumbbellPlot(before_campaign, after_campaign, ...
        'MarkerSize', importance, ...           
        'labelX1', 'Pre-Campaign', ...
        'labelX2', 'Post-Campaign', ...
        'Title', 'Marketing Campaign Impact (marker size = revenue weight)', ...
        'YLabels', products);
%}

%{
    %% Example of plot with background bands and text inside the dots.
    X1 = [45, 62, 38, 71, 55, 48];
    X2 = [58, 55, 75, 84, 68, 52];
    categories = {'Sales', 'Marketing', 'IT', 'HR', 'Finance', 'Operations'};

    figure;
    dumbbellPlot(X1,X2, "YLabels", categories, "Title", "Horizontal Plot with background bands", ...
        "TextInside", true, "labelX1","Q1", "labelX2", "Q2", "Background", "bands");
%}

%% Beginning of code
% import the object from +graphics/DumbellChart.m
import graphics.DumbbellChart

plotType= "single"; %default
plotTypeidx= find(strcmp(varargin, "plotType"));
if ~isempty(plotTypeidx)
    if plotTypeidx + 1 <= length(varargin)
        plotType = varargin{plotTypeidx + 1};
    else
        error("dumbbellPlot:MissingPlotTypeValue", "plotType 'name',value argument requires a value")
    end
end
%% Required input validation
if nargin < 1
    error("dumbbellPlot:MissingRequiredInput", "Missing required input argument X")
end

% check X, it can be vector, or table
if istable(X)
    switch width(X)
        case 1
            if length(varargin) < 1
                error("dumbbellPlot:MissingRequiredInput", "When X is a single-column table, a second set of values must be provided")
            end

            X1= X{:,1};
            labelX1= X.Properties.VariableNames{1};

            if istable(varargin{1})
                if width(varargin{1}) ~= 1
                    error("dumbbellPlot:InvalidTableWidth", ...
                        "Second table argument must have exactly one column")
                end
                X2 = varargin{1}{:,1};
                labelX2= varargin{1}.Properties.VariableNames{1};
                vararginStart = 2;

            elseif isnumeric(varargin{1}) && isvector(varargin{1})
                if length(varargin{1}) ~= length(X1)
                    error("dumbbellPlot:DimensionMismatch", "X2 length must match X1 length")
                end
                X2 = varargin{1};
                labelX2= "X2";
                vararginStart = 2;
            else
                error("dumbbellPlot:WrongInput", ...
                    "X2 must be a single-column table or numeric vector with compatible size")
            end

        case 2
            if plotType ~= "single"
                error("dumbbellPlot:WrongInput", "Tables with 2 variables are only accepted if plotType is 'single'")
            end

            X1 = X{:,1};
            X2 = X{:,2};
            labelX1= X.Properties.VariableNames{1};
            labelX2= X.Properties.VariableNames{2};
            vararginStart= 1;

        otherwise %3 or more columns
            if plotType == "double"
                if width(X) == 4
                    X1 = X{:,1};
                    X2 = X{:,2};
                    X3 = X{:,3};
                    X4 = X{:,4};

                    % not sure if to leave default, TO BE DECIDED
                    labelX1 = X.Properties.VariableNames{1};
                    labelX2 = X.Properties.VariableNames{2};
                    vararginStart = 1;
                else
                    if length(varargin) < 4
                        error("dumbbellPlot:MissingTableVariables", ...
                            "When X is a multi-column table, you must specify 4 variable names for Value1, 2, 3 and 4 in a double plotType")
                    end



                    varNames = {varargin{1}, varargin{2}, varargin{3}, varargin{4}};
                    for i = 1:4
                        if ~(ischar(varNames{i}) || isstring(varNames{i}))
                            error("dumbbellPlot:InvalidVariableName", ...
                                "Variable name %d must be a string or character array", i);
                        end
                        if ~ismember(varNames{i}, X.Properties.VariableNames)
                            error("dumbbellPlot:VariableNotFound", ...
                                "Variable '%s' not found in table", varNames{i});
                        end
                    end

                    vararginStart= 5;

                    % assign data set values
                    X1= X.(varNames{1});
                    X2= X.(varNames{2});
                    X3= X.(varNames{3});
                    X4= X.(varNames{4});

                    % not sure if leaving default is better DECIDE LATER
                    labelX1 = varNames{1};
                    labelX2 = varNames{2};

                end

            elseif plotType == "single"
                if length(varargin) < 2
                    error("dumbbellPlot:MissingTableVariables", ...
                        "When X is a multi-column table, you must specify two variable names for Value1 and Value2")
                end

                varName1 = varargin{1};
                varName2 = varargin{2};
                vararginStart= 3;

                if ~(ischar(varName1) || isstring(varName1))
                    error("dumbbellPlot:InvalidVariableName", "Variable name 1 must be a string or character array")
                end
                if ~(ischar(varName2) || isstring(varName2))
                    error("dumbbellPlot:InvalidVariableName", "Variable name 2 must be a string or character array")
                end

                if ~ismember(varName1, X.Properties.VariableNames)
                    error("dumbbellPlot:VariableNotFound", "Variable '%s' not found in table", varName1)
                end
                if ~ismember(varName2, X.Properties.VariableNames)
                    error("dumbbellPlot:VariableNotFound", "Variable '%s' not found in table", varName2)
                end

                X1= X.(varName1);
                X2= X.(varName2);
                labelX1 = varName1;
                labelX2 = varName2;

            end
    end

    % add Yticklabels from table or set default ones
    if ~isempty(X.Properties.RowNames)
        YLabels = X.Properties.RowNames;
    else
        YLabels = cellstr("Row " + string(1:length(X1)));
    end

elseif isvector(X) && isnumeric(X)
    if length(varargin) < 1
        error("dumbbellPlot:MissingRequiredInput", "If X1 is a vector, a second vector is required")
    end
    if ~(isvector(varargin{1}) && isnumeric(varargin{1}))
        error("dumbbellPlot:MissingRequiredInput", "Second input argument is not a numeric vector")
    end

    X1= X;
    X2= varargin{1};
    vararginStart= 2;
    if length(X1) ~= length(X2)
        error("dumbbellPlot:DimensionMismatch", ...
            "All data vectors must have the same length");
    end

    labelX1 = "X1";
    labelX2 = "X2";
    YLabels = cellstr("Row" + string(1:length(X1)));
    if plotType == "double"
        if length(varargin) < 3
            error("dumbbellPlot:MissingRequiredInput", "If plotType is 'double', 4 sets of data must be inserted")
        end
        if ~(isvector(varargin{2}) && isnumeric(varargin{2}))
            error("dumbbellPlot:MissingRequiredInput", "Third input argument is not a numeric vector")
        end
        if ~(isvector(varargin{3}) && isnumeric(varargin{3}))
            error("dumbbellPlot:MissingRequiredInput", "Fourth input argument is not a numeric vector")
        end

        X3= varargin{2};
        X4= varargin{3};
        vararginStart= 4;

        if length(X3) ~= length(X1) || length(X4) ~= length(X1)
            error("dumbbellPlot:DimensionMismatch", ...
                "All data vectors must have the same length");
        end
    end
else
    error("dumbbellPlot:InvalidInputType", "X must be a table or a numeric vector")
end

%% optional inputs

% struct required by publishFS function to create the help file
% default values with a decimal .1 are used as sentinels to detect if
% default value was used.
options= struct("labelX1", labelX1, ...
                "labelX2", labelX2, ...
                "plotType", plotType, ... 
                "Title", "", ...
                "YLabels", {YLabels}, ... 
                "orientation", "horizontal", ...
                "Color", getColorPalette("default"), ...
                "MarkerSize", 150.1, ...
                "LineWidth", 3.1, ...
                "TextSize", 12.1, ...
                "TextInside", false, ...
                "ColorDist", "false", ...
                "AxesFontSize", 13, ...
                "Background", "none");

UserOptions = struct();
if length(varargin) >= vararginStart
    for i = vararginStart:2:length(varargin)
        if i+1 <= length(varargin)
            UserOptions.(varargin{i}) = varargin{i+1};
        end
    end
end

if isfield(UserOptions, "labelX1")
    validateattributes(UserOptions.labelX1, {'string', 'char'}, {}, 'dumbellPlot', 'labelX1')
    options.labelX1 = string(UserOptions.labelX1);
end

if isfield(UserOptions, "labelX2")
    validateattributes(UserOptions.labelX2, {'string', 'char'}, {}, 'dumbellPlot', 'labelX2')
    options.labelX2 = string(UserOptions.labelX2);
end

if options.plotType == "double"
    options.Title = ["Year 1"; "Year 2"]; % set default in case of double plot
end

if isfield(UserOptions, "Title")
    validateattributes(UserOptions.Title, {'string', 'char', 'cell'}, {}, 'dumbellPlot', 'Title')

    if iscell(UserOptions.Title)
        tempTitle = UserOptions.Title;
    else
        tempTitle = cellstr(UserOptions.Title);
    end

    if options.plotType == "single"
        if numel(tempTitle) ~= 1
            warning("dumbbellPlot:TitleMismatch", "Single plot requires 1 title, default will be used")
        else
            options.Title=string(tempTitle{1});
        end
    elseif options.plotType == "double"
        if numel(tempTitle) ~= 2
            warning("dumbbellPlot:TitleMismatch", "Double plot requires 2 titles, default will be used")
        else
            options.Title=tempTitle;
        end
    end
end

if isfield(UserOptions, "YLabels")
    validateattributes(UserOptions.YLabels, {'string', 'char', 'cell'}, {}, 'dumbbellPlot', 'YLabels')

    if iscell(UserOptions.YLabels)
        tempYLabels = UserOptions.YLabels;
    else
        tempYLabels = cellstr(UserOptions.YLabels);
    end

    if length(tempYLabels) ~= length(X1)
        warning("dumbbellPlot:YLabelsMismatch", "Ylabels and X lenght does not match, default YLabels will be used")
    else
        options.YLabels= tempYLabels;
    end
end

if isfield(UserOptions, "orientation")
    validateattributes(UserOptions.orientation, {'string', 'char'}, {}, 'dumbbellPlot', 'orientation')
    options.orientation= UserOptions.orientation;
end

if isfield(UserOptions, "Color")
    if (ischar(UserOptions.Color) && isrow(UserOptions.Color)) || (isstring(UserOptions.Color) && isscalar(UserOptions.Color))
        options.Color = getColorPalette(UserOptions.Color);
    else
        try
            options.Color = validatecolor(UserOptions.Color, "multiple");
            if size(options.Color, 1) ~= 2
                warning("dumbbellPlot:InvalidColor", "Color must provide exactly 2 valid colors, default palette will be used")
                options.Color = getColorPalette("default");
            end
        catch ME
            warning("dumbbellPlot:InvalidColor", "Invalid Color, default palette will be used")
            options.Color = getColorPalette("default");
        end
    end
end

if isfield(UserOptions, "MarkerSize")
    if isnumeric(UserOptions.MarkerSize) && isscalar(UserOptions.MarkerSize)
        options.MarkerSize= UserOptions.MarkerSize;
    elseif isnumeric(UserOptions.MarkerSize) && length(UserOptions.MarkerSize) == length(X1)
        options.MarkerSize= UserOptions.MarkerSize;
    else
        warning("dumbbellPlot:InvalidMarkerSize", "Invalid Marker size provided, using default")
    end
end

if isfield(UserOptions, "LineWidth")
    if isnumeric(UserOptions.LineWidth) && isscalar(UserOptions.LineWidth)
        options.LineWidth= UserOptions.LineWidth;
    else
        warning("dumbbellPlot:InvalidLineWidth", "Invalid Line Width provided, using default")
        options.LineWidth= 3.1;
    end
end

if isfield(UserOptions, "TextSize")
    if isnumeric(UserOptions.TextSize) && isscalar(UserOptions.TextSize)
        options.TextSize= UserOptions.TextSize;
    else
        warning("dumbbellPlot:InvalidTextSize", "Invalid Text size provided, using default")
    end
end

if isfield(UserOptions, "TextInside")
    if islogical(UserOptions.TextInside) && isscalar(UserOptions.TextInside)
        options.TextInside= UserOptions.TextInside;
    else
        warning("dumbbellPlot:InvalidTextInside", "Text Inside needs to be a scalar logical (true or false) value!")
    end
end

if isfield(UserOptions, "ColorDist")
    if (ischar(UserOptions.ColorDist) && isrow(UserOptions.ColorDist)) || (isstring(UserOptions.ColorDist) && isscalar(UserOptions.ColorDist))
        options.ColorDist= string(UserOptions.ColorDist);

        % validate type
        validTypes = ["false", "directional", "magnitude", "robust"];
        if ~ismember(options.ColorDist, validTypes)
            warning("dumbbellPlot:InvalidColorDist", ...
                "ColorDist must be 'false', 'directional', 'magnitude', or 'robust'. Using default (false)")   
        end
    else
        warning("dumbbellPlot:InvalidColorDist", "ColorDist needs to be char or string")
    end
end

if isfield(UserOptions, "AxesFontSize")
    if isnumeric(UserOptions.AxesFontSize) && isscalar(UserOptions.AxesFontSize)
        options.AxesFontSize= UserOptions.AxesFontSize;
    else
        warning("dumbbellPlot:InvalidFontSize","Invalid Font Size for axes provided, using default")
    end
end

if isfield(UserOptions, "Background")
    validBG= ["bands", "none"];
    if ismember(UserOptions.Background, validBG)
        options.Background= string(UserOptions.Background);
    else
        warning("dumbbellPlot:InvalidBakgroundOption","Invalid Background option provided, using default (no background)")
    end
else
end
%% main body
switch strcat(options.plotType,"_",options.orientation)
    case "single_horizontal"
        
        ax = gca;
        chart = DumbbellChart(X1,X2,options.YLabels,options.Color,options.MarkerSize, ...
                            options.LineWidth,options.TextSize,options.TextInside, ...
                            options.ColorDist,options.AxesFontSize,options.Background);
        [h1, h2] = chart.build(ax);
        legend([h1 h2], {options.labelX1, options.labelX2}, "Location","best")

        if options.ColorDist ~= "false"
            colormap(ax, turbo(256))
            cb = colorbar(ax);
            cb.Label.String = "Difference: "+options.ColorDist;

            if options.ColorDist == "directional"
                clim(ax, [-1 1]); % relative scale
                cb.Label.String = "Difference: Directional (relative scale)";
            elseif options.ColorDist == "robust"
                clim(ax, [-2, 2]); % remember to match value that used in clipping
            end
        end

        % title
        if options.Title ~= ""
            title(ax, options.Title);
        end

    case "single_vertical"
        ax = gca;
        chart = DumbbellChart(X1,X2,options.YLabels,options.Color,options.MarkerSize, ...
                            options.LineWidth,options.TextSize,options.TextInside, ...
                            options.ColorDist,options.AxesFontSize,options.Background);
        [h1, h2] = chart.buildVertical(ax);
        legend([h1 h2], {options.labelX1, options.labelX2}, "Location","best")

        if options.ColorDist ~= "false"
            colormap(ax, turbo(256))
            cb = colorbar(ax);
            cb.Label.String = "Difference: "+options.ColorDist;

            if options.ColorDist == "directional"
                clim(ax, [-1 1]); 
                cb.Label.String = "Difference: Directional (relative scale)";
            elseif options.ColorDist == "robust"
                clim(ax, [-2, 2]); % remember to match clipping
            end
        end

        % title
        if options.Title ~= ""
            title(ax, options.Title);
        end

    case "double_horizontal"
        tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        ax1 = nexttile;
        chart = DumbbellChart(X1,X2,options.YLabels,options.Color,options.MarkerSize, ...
                            options.LineWidth,options.TextSize,options.TextInside, ...
                            options.ColorDist,options.AxesFontSize,options.Background);
        [h1, h2] = chart.build(ax1);
        legend(ax1, [h1 h2], {options.labelX1, options.labelX2}, "Location","best")
        
        ax2= nexttile;
        chart2 = DumbbellChart(X3,X4,options.YLabels,options.Color,options.MarkerSize, ...
                            options.LineWidth,options.TextSize,options.TextInside, ...
                            options.ColorDist,options.AxesFontSize,options.Background);
        [~, ~] = chart2.build(ax2);

        title(ax1, options.Title{1});
        title(ax2, options.Title{2});

        if options.ColorDist ~= "false"
            colormap(ax1, turbo(256))
            colormap(ax2, turbo(256))
            cb = colorbar(ax2);
            cb.Layout.Tile = 'east';
            cb.Label.String = "Difference: "+options.ColorDist;

            if options.ColorDist == "directional"
                clim(ax2, [-1 1]);
                cb.Label.String = "Difference: Directional (relative scale)";
            elseif options.ColorDist == "robust"
                clim(ax2, [-2, 2]); % remember to match clipping
            end
        end

        % return value
        ax = [ax1; ax2];

    case "double_vertical"
        tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        ax1 = nexttile;
        chart = DumbbellChart(X1,X2,options.YLabels,options.Color,options.MarkerSize, ...
                            options.LineWidth,options.TextSize,options.TextInside, ...
                            options.ColorDist,options.AxesFontSize,options.Background);
        [h1, h2] = chart.buildVertical(ax1);
        
        ax2 = nexttile;
        chart2 = DumbbellChart(X3,X4,options.YLabels,options.Color,options.MarkerSize, ...
                            options.LineWidth,options.TextSize,options.TextInside, ...
                            options.ColorDist,options.AxesFontSize,options.Background);
        [~, ~] = chart2.buildVertical(ax2);
        legend(ax2, [h1 h2], {options.labelX1, options.labelX2}, "Location","best")
        ax2.YAxisLocation= "right";
        ax2.YLabel.String= "";
        
        title(ax1, options.Title{1});
        title(ax2, options.Title{2});

        if options.ColorDist ~= "false"
            colormap(ax1, turbo(256))
            colormap(ax2, turbo(256))
            cb = colorbar(ax2);
            cb.Layout.Tile = 'east';
            cb.Label.String = "Difference: "+options.ColorDist;

            if options.ColorDist == "directional"
                clim(ax2, [-1 1]);
                cb.Label.String = "Difference: Directional (relative scale)";
            elseif options.ColorDist == "robust"
                clim(ax2, [-2, 2]); % remember to match clipping value
            end
        end

        % return value
        ax = [ax1; ax2];
end
    
%% helper function to support pre-made color palettes
    function colors = getColorPalette(paletteName)
        switch paletteName
            case "default"
                colors = [0.92, 0.78, 0.38; % Sand yellow
                    0.17, 0.44, 0.26]; %Forest green
            case "colorblind"
                colors = [0.90, 0.62, 0.00; % Amber orange
                    0.00, 0.45, 0.70]; % Ocean blue
            case "ruby_jade"
                colors = [0.78, 0.15, 0.28; % Ruby red
                    0.25, 0.60, 0.50]; % Jade green
            case "cherry_sky"
                colors = [0.85, 0.25, 0.38; % Cherry red
                    0.52, 0.75, 0.92]; % Sky blue
            case "red_blue"
                colors = [0.88, 0.30, 0.30; % Coral Red
                    0.20, 0.42, 0.85]; % Royal blue
            otherwise
                warning("dumbbellPlot:UnknownPalette", ...
                    "Unknown palette '%s', using default", paletteName);
                colors = [0.92, 0.78, 0.38; 0.17, 0.44, 0.26];
        end
    end
end
%FScategory:VIS-Mult