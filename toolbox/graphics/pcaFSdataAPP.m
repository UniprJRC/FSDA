classdef pcaFSdataAPP < matlab.apps.AppBase
    %pcaFSdataAPP Load a data file and launch pcaFS.

    properties (Access = public)
        PCAFSDataImportUIFigure  matlab.ui.Figure
        GridLayout               matlab.ui.container.GridLayout
        ControlPanel             matlab.ui.container.Panel
        ControlGrid              matlab.ui.container.GridLayout
        LoadButton               matlab.ui.control.Button
        RunButton                matlab.ui.control.Button
        FileLabel                matlab.ui.control.Label
        DataSizeLabel            matlab.ui.control.Label
        NumericColumnsLabel      matlab.ui.control.Label
        ForwardedOptionsLabel    matlab.ui.control.Label
        StatusLabel              matlab.ui.control.Label
        PreviewTable             matlab.ui.control.Table
    end

    properties (Access = private)
        DataTable table = table()
        FileName char = ''
        ForwardedOptions cell = {}
    end

    methods (Access = private)
        function startupFcn(app, varargin)
            [initialFile, app.ForwardedOptions] = parseStartupInputs(app, varargin);
            updateForwardedOptionsLabel(app)
            updateStatus(app, 'Choose a data file to begin.', false)

            if ~isempty(initialFile)
                loadFile(app, initialFile)
            end
        end

        function [initialFile, forwardedOptions] = parseStartupInputs(app, inputs)
            initialFile = '';
            forwardedOptions = inputs;

            if numel(inputs) >= 2 && isTextScalar(app, inputs{1}) && ...
                    any(strcmpi(char(inputs{1}), {'InitialFile', 'FileName'}))
                initialFile = char(inputs{2});
                forwardedOptions = inputs(3:end);
            elseif isscalar(inputs) && isTextScalar(app, inputs{1}) && isfile(char(inputs{1}))
                initialFile = char(inputs{1});
                forwardedOptions = {};
            end

            if ~isempty(forwardedOptions)
                [forwardedOptions{:}] = convertStringsToChars(forwardedOptions{:});
            end
        end

        function LoadButtonPushed(app, ~)
            filters = { ...
                '*.xlsx;*.xls;*.xlsm;*.csv;*.txt;*.dat', 'Data files (*.xlsx, *.xls, *.csv, *.txt, *.dat)'; ...
                '*.xlsx;*.xls;*.xlsm', 'Excel files (*.xlsx, *.xls, *.xlsm)'; ...
                '*.csv', 'CSV files (*.csv)'; ...
                '*.txt;*.dat', 'Text files (*.txt, *.dat)'; ...
                '*.*', 'All files (*.*)'};
            [fileName, pathName] = uigetfile(filters, 'Select data file for pcaFS');
            if isequal(fileName, 0)
                updateStatus(app, 'File selection cancelled.', false)
                return
            end
            loadFile(app, fullfile(pathName, fileName))
        end

        function RunButtonPushed(app, ~)
            if isempty(app.DataTable) || width(app.DataTable) == 0
                uialert(app.PCAFSDataImportUIFigure, 'Load a data file before running pcaFS.', 'No data loaded')
                return
            end

            if mod(numel(app.ForwardedOptions), 2) ~= 0
                uialert(app.PCAFSDataImportUIFigure, ...
                    'Optional arguments must be supplied as name-value pairs.', ...
                    'Invalid varargin')
                return
            end

            app.RunButton.Enable = 'off';
            updateStatus(app, 'Calling pcaFS...', false)
            drawnow

            try
                assignin('base', 'pcaFSdataAPP_loadedData', app.DataTable)
                assignin('base', 'pcaFSdataAPP_varargin', app.ForwardedOptions)
                out = pcaFS(app.DataTable, app.ForwardedOptions{:});
                assignin('base', 'pcaFSdataAPP_output', out)
                updateStatus(app, 'pcaFS completed.', false)
            catch ME
                updateStatus(app, ME.message, true)
                uialert(app.PCAFSDataImportUIFigure, ME.message, 'pcaFS failed')
            end

            app.RunButton.Enable = 'on';
        end

        function loadFile(app, fileName)
            try
                fileName = char(fileName);
                T = readDataFile(app, fileName);
                [Tnumeric, detail] = prepareDataForPCA(app, T);

                app.DataTable = Tnumeric;
                app.FileName = fileName;
                [~, name, ext] = fileparts(fileName);
                app.FileLabel.Text = [name ext];
                app.FileLabel.Tooltip = {fileName};
                app.DataSizeLabel.Text = sprintf('%d rows, %d numeric variables', height(Tnumeric), width(Tnumeric));
                app.NumericColumnsLabel.Text = strjoin(Tnumeric.Properties.VariableNames, ', ');
                app.NumericColumnsLabel.Tooltip = {app.NumericColumnsLabel.Text};
                app.PreviewTable.Data = Tnumeric;
                app.RunButton.Enable = 'on';

                if strlength(detail) > 0
                    updateStatus(app, char(detail), false)
                else
                    updateStatus(app, 'Data loaded. Ready to call pcaFS.', false)
                end
            catch ME
                app.DataTable = table();
                app.FileName = '';
                app.FileLabel.Text = 'No file selected';
                app.FileLabel.Tooltip = {};
                app.DataSizeLabel.Text = 'No data loaded';
                app.NumericColumnsLabel.Text = '';
                app.NumericColumnsLabel.Tooltip = {};
                app.PreviewTable.Data = table();
                app.RunButton.Enable = 'off';
                updateStatus(app, ME.message, true)
                uialert(app.PCAFSDataImportUIFigure, ME.message, 'Import failed')
            end
        end

        function T = readDataFile(~, fileName)
            if ~isfile(fileName)
                error('FSDA:pcaFSdataAPP:FileNotFound', 'File not found: %s', fileName)
            end

            [~, ~, ext] = fileparts(fileName);
            ext = lower(ext);
            try
                switch ext
                    case {'.xls', '.xlsx', '.xlsm', '.xlsb'}
                        opts = detectImportOptions(fileName, 'FileType', 'spreadsheet');
                    otherwise
                        opts = detectImportOptions(fileName, 'FileType', 'text');
                end
                T = readtable(fileName, opts);
            catch
                T = readtable(fileName);
            end

            if height(T) == 0 || width(T) == 0
                error('FSDA:pcaFSdataAPP:EmptyData', 'The selected file did not contain a non-empty table.')
            end
        end

        function [Tnumeric, detail] = prepareDataForPCA(app, T)
            isNumeric = false(1, width(T));
            for j = 1:width(T)
                values = T{:, j};
                isNumeric(j) = (isnumeric(values) || islogical(values)) && ismatrix(values) && size(values, 2) == 1;
            end

            if ~any(isNumeric)
                error('FSDA:pcaFSdataAPP:NoNumericColumns', 'The selected file does not contain numeric columns for pcaFS.')
            end

            Tnumeric = T(:, isNumeric);
            for j = 1:width(Tnumeric)
                varName = Tnumeric.Properties.VariableNames{j};
                Tnumeric.(varName) = double(Tnumeric{:, j});
            end

            rowNames = rowNamesFromFirstTextColumn(app, T, isNumeric);
            if ~isempty(rowNames)
                Tnumeric.Properties.RowNames = rowNames;
            end

            X = table2array(Tnumeric);
            keepRows = all(isfinite(X), 2);
            removedRows = sum(~keepRows);
            if removedRows > 0
                Tnumeric = Tnumeric(keepRows, :);
            end

            if height(Tnumeric) < 2
                error('FSDA:pcaFSdataAPP:TooFewRows', 'pcaFS requires at least two valid observations.')
            end
            if width(Tnumeric) < 2
                error('FSDA:pcaFSdataAPP:TooFewColumns', 'pcaFS requires at least two numeric variables for the dynamic biplot.')
            end

            skippedColumns = sum(~isNumeric);
            parts = strings(1, 0);
            if skippedColumns > 0
                parts(end + 1) = sprintf('%d nonnumeric column(s) skipped.', skippedColumns);
            end
            if removedRows > 0
                parts(end + 1) = sprintf('%d row(s) with missing or infinite values removed.', removedRows);
            end

            if isempty(parts)
                detail = "";
            else
                detail = strjoin(parts, ' ');
            end
        end

        function rowNames = rowNamesFromFirstTextColumn(~, T, isNumeric)
            rowNames = {};
            for j = find(~isNumeric)
                values = T{:, j};
                if isstring(values) || iscellstr(values) || iscategorical(values) || isdatetime(values) || isduration(values)
                    names = string(values);
                elseif iscell(values)
                    try
                        names = string(values);
                    catch
                        continue
                    end
                else
                    continue
                end

                if any(ismissing(names))
                    continue
                end
                names = strtrim(names);
                if all(strlength(names) > 0) && numel(unique(names)) == height(T)
                    rowNames = cellstr(names);
                    return
                end
            end
        end

        function updateForwardedOptionsLabel(app)
            if isempty(app.ForwardedOptions)
                app.ForwardedOptionsLabel.Text = 'No pcaFS optional arguments supplied';
            else
                app.ForwardedOptionsLabel.Text = ['pcaFS varargin: ' formatVarargin(app, app.ForwardedOptions)];
            end
            app.ForwardedOptionsLabel.Tooltip = {app.ForwardedOptionsLabel.Text};
        end

        function txt = formatVarargin(~, args)
            txtParts = strings(1, numel(args));
            for k = 1:numel(args)
                value = args{k};
                if ischar(value)
                    txtParts(k) = "'" + string(value) + "'";
                elseif isstring(value) && isscalar(value)
                    txtParts(k) = """" + value + """";
                elseif isnumeric(value) || islogical(value)
                    if isscalar(value)
                        txtParts(k) = string(value);
                    else
                        txtParts(k) = sprintf('%s array', class(value));
                    end
                else
                    txtParts(k) = sprintf('%s value', class(value));
                end
            end
            txt = char(strjoin(txtParts, ', '));
        end

        function updateStatus(app, message, isError)
            app.StatusLabel.Text = message;
            if isError
                app.StatusLabel.FontColor = [0.65 0 0];
            else
                app.StatusLabel.FontColor = [0.15 0.15 0.15];
            end
            drawnow limitrate
        end

        function tf = isTextScalar(~, value)
            tf = (ischar(value) && (isrow(value) || isempty(value))) || ...
                (isstring(value) && isscalar(value));
        end

        function createComponents(app)
            app.PCAFSDataImportUIFigure = uifigure('Visible', 'off');
            app.PCAFSDataImportUIFigure.Position = [100 100 860 500];
            app.PCAFSDataImportUIFigure.Name = 'pcaFS data import';

            app.GridLayout = uigridlayout(app.PCAFSDataImportUIFigure, [1 2]);
            app.GridLayout.ColumnWidth = {280, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.Padding = [12 12 12 12];
            app.GridLayout.ColumnSpacing = 12;

            app.ControlPanel = uipanel(app.GridLayout);
            app.ControlPanel.Title = 'Data file';
            app.ControlPanel.Layout.Row = 1;
            app.ControlPanel.Layout.Column = 1;

            app.ControlGrid = uigridlayout(app.ControlPanel, [9 1]);
            app.ControlGrid.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', '1x', 'fit', 'fit', 'fit'};
            app.ControlGrid.ColumnWidth = {'1x'};
            app.ControlGrid.Padding = [10 10 10 10];
            app.ControlGrid.RowSpacing = 8;

            app.LoadButton = uibutton(app.ControlGrid, 'push');
            app.LoadButton.Text = 'Load data file';
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Layout.Row = 1;
            app.LoadButton.Layout.Column = 1;

            app.FileLabel = uilabel(app.ControlGrid);
            app.FileLabel.Text = 'No file selected';
            app.FileLabel.WordWrap = 'on';
            app.FileLabel.Layout.Row = 2;
            app.FileLabel.Layout.Column = 1;

            app.DataSizeLabel = uilabel(app.ControlGrid);
            app.DataSizeLabel.Text = 'No data loaded';
            app.DataSizeLabel.Layout.Row = 3;
            app.DataSizeLabel.Layout.Column = 1;

            app.NumericColumnsLabel = uilabel(app.ControlGrid);
            app.NumericColumnsLabel.Text = '';
            app.NumericColumnsLabel.WordWrap = 'on';
            app.NumericColumnsLabel.Layout.Row = 4;
            app.NumericColumnsLabel.Layout.Column = 1;

            app.ForwardedOptionsLabel = uilabel(app.ControlGrid);
            app.ForwardedOptionsLabel.Text = '';
            app.ForwardedOptionsLabel.WordWrap = 'on';
            app.ForwardedOptionsLabel.Layout.Row = 5;
            app.ForwardedOptionsLabel.Layout.Column = 1;

            app.RunButton = uibutton(app.ControlGrid, 'push');
            app.RunButton.Text = 'Run pcaFS';
            app.RunButton.Enable = 'off';
            app.RunButton.ButtonPushedFcn = createCallbackFcn(app, @RunButtonPushed, true);
            app.RunButton.Layout.Row = 8;
            app.RunButton.Layout.Column = 1;

            app.StatusLabel = uilabel(app.ControlGrid);
            app.StatusLabel.Text = '';
            app.StatusLabel.WordWrap = 'on';
            app.StatusLabel.Layout.Row = 9;
            app.StatusLabel.Layout.Column = 1;

            app.PreviewTable = uitable(app.GridLayout);
            app.PreviewTable.Layout.Row = 1;
            app.PreviewTable.Layout.Column = 2;

            app.PCAFSDataImportUIFigure.Visible = 'on';
        end
    end

    methods (Access = public)
        function app = pcaFSdataAPP(varargin)
            createComponents(app)
            registerApp(app, app.PCAFSDataImportUIFigure)
            runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            if nargout == 0
                clear app
            end
        end

        function delete(app)
            if ~isempty(app.PCAFSDataImportUIFigure) && isvalid(app.PCAFSDataImportUIFigure)
                delete(app.PCAFSDataImportUIFigure)
            end
        end
    end
end
