function [latex_string , disp_string , latex_string_full] = tabledisp(T, precision, stylerc, filename)
%tabledisp displays in good format a table or array in the command line
%
%<a href="matlab: docsearchFS('tabledisp')">Link to the help function</a>
%
%  tabledisp displays a table or o bi-dimensional array in the MATLAB
%  command line or in a uitable and also in an excel file.
%
%  Required input arguments:
%
%            T: Input table. Table or bi-dimensional array. The table/array
%               contains a data matrix with n observations on v variables.
%               Missing values (NaN's) and infinite values (Inf's) are allowed.
%
%  Optional input arguments:
%
%    precision: Number of digits after the comma, to display.
%               Scalar. An integer value.
%               Example - 'precision',2
%               Data Types - single|double
%
%    stylerc:   The display style of the raw and column names. Character
%               array. It can take these values:
%               'b' which stands for bold;
%               'i' which stands for italic;
%               'n' which stand for normal (the default).
%               If the user mis-specify the style, 'n is used.
%               Example - 'stylerc','b'
%               Data Types - char
%
%    filename:  The name of a file where the formatted table will be saved.
%               Array of characters. A filename without empty spaces. Note
%               that the latex string in the file contains the proper \n
%               characters and can be easily incorporated in latex
%               documents. If filename is 'excel', then an excel file is
%               generated with the fixed name 'tabledisp_excel.xlsx'.
%               Example - 'filename','pippo.txt'
%               Data Types - char
%  Output:
%
%         latex_string :   Latex string with table display. Character array.
%                          The string can be used in uitables, annotation
%                          statements or other matlab tools receving latex
%                          sentences. Note that this string by default has
%                          no \n characters, as they are not treated by the
%                          annotation statement. Note also that this string
%                          is automatically and appropriately cut to comply
%                          with the 1200 characters limitaion of
%                          annotations and latex statements in graphical
%                          objects.
%          disp_string :   String with table display. Character array.
%                          The string can be displayed in the command
%                          window with fprintf(disp_string).
%    latex_string_full :   Latex string with table display. Character array.
%                          This is the full version of latex_string,
%                          without the cut on 1200 characters. It is
%                          reported in output in case there is a need to
%                          copy-paste the full table in LaTeX documents. 
%
% See also: disp
%
% References:
%
%    Lamport, L. (1994), LATEX: a document preparation system:
%    user's guide and reference manual. Addison-Wesley Longman Publishing
%    Co., Inc., USA.
%
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('tabledisp')">Link to the help function</a>
%
%$LastChangedDate:: 2018-01-29 18:52:14 #$: Date of the last commit

% Examples:

%{
    %% Dispay a table in MATLAB annotations.
    % load a generic dataset
    load patients;
    % put it in a table
    T = table(Age,Height,Weight,Systolic,Diastolic, ...
              'RowNames',LastName);
    % keep just few lines, as example
    Tsmall = T(1:10,:);
    % run tabledisp
    [latex_string , disp_string] = tabledisp(Tsmall);

    % Now use the output of tabledisp to diplay the latex table in the
    % annotation of a figure
    hf = figure;
    annotation(hf,'Textbox','String',latex_string,...
        'FitBoxToText','on','Interpreter','latex',...
        'FontName',get(0,'FixedWidthFontName'),'FontSize',14,...
        'Units','Normalized','Position',[0 0 1 1]);

    % An this is to display the table in the command window without disp
    fprintf(disp_string);
%}

%{
    %% Dispay an array in a MATLAB annotation, with 6 digits of precision.
    % generate data
    X = 1000*randn(15,5);
    % and run tabledisp on them with the specification of the precision
    [latex_string , disp_string] = tabledisp(X,6);

    % This is to produce a latex table display in a figure
    hf = figure;
    annotation(hf,'Textbox','String',latex_string,...
        'FitBoxToText','on','Interpreter','latex',...
        'FontName',get(0,'FixedWidthFontName'),'FontSize',14,...
        'Units','Normalized','Position',[0 0 1 1]);

    % This is to display the table in the command window without disp
    fprintf(disp_string);
%}

%{
    %% Dispay an array in a MATLAB annotation with bold row/column names.
    close all
    % generate data
    X = 1000*randn(10,5);
    % and run tabledisp on them with the specification of the italic style
    [latex_string_i , disp_string] = tabledisp(X,[],'i');

    hfi = figure;
    annotation(hfi,'Textbox','String',latex_string_i,...
        'FitBoxToText','on','Interpreter','latex',...
        'FontName',get(0,'FixedWidthFontName'),'FontSize',14,...
        'Units','Normalized','Position',[0 0 1 1]);

    % and now bold style
    [latex_string_b , disp_string] = tabledisp(X,[],'b');

    hfb = figure;
    annotation(hfb,'Textbox','String',latex_string_b,...
        'FitBoxToText','on','Interpreter','latex',...
        'FontName',get(0,'FixedWidthFontName'),'FontSize',14,...
        'Units','Normalized','Position',[0 0 1 1]);

    cascade;
%}

%{
    %% Dispay a large table in a MATLAB annotation.
    % generate data
    X = 1000*randn(150,5);
    % and run tabledisp on them with the specification of the precision
    [latex_string , disp_string , latex_string_full] = tabledisp(X,6,'b');

    % This is to produce a latex table display in a figure
    hf = figure;
    annotation(hf,'Textbox','String',latex_string,...
        'FitBoxToText','on','Interpreter','latex',...
        'FontName',get(0,'FixedWidthFontName'),'FontSize',14,...
        'Units','Normalized','Position',[0 0 1 1]);

    % This is to display the table in the command window without disp
    fprintf(disp_string);
%}

%{
    % Dispay an array in a MATLAB annotation, and save it in a specific file.
    X = randn(10,5);
    % and run tabledisp on them with the specification of the precision
    [latex_string , disp_string] = tabledisp(X,2,'','myfile.txt');

%}

%{
    % Save an array in a excel file.
    X = randn(10,5);
    % and run tabledisp on them with the specification of the precision
    % An excel file named tabledisp_excel.xlsx is created in the current folder.
    [latex_string , disp_string] = tabledisp(X,2,'','excel');

%}

%% Beginning of code

% DEFAULT OPTIONAL PARAMETERS

% filename to save output
if nargin == 4 && strcmp(filename,'excel')
    filename  = ['tabledisp_' filename];
end
if nargin < 4 || ~ischar(filename)
    filename  = '';
end

% style of the row/column names
if nargin < 3 || isempty(stylerc) || ~ismember(stylerc , {'n','i','b'})
    stylerc = 'n';
end
switch stylerc
    case 'b'
        stylerclatex  = '\\bf ';  
        stylerclatex2 = '\bf ';   % note: \bf vs \\bf
    case 'i'
        stylerclatex  = '\\it ';
        stylerclatex2 = '\it ';   % note: \it vs \\it
    otherwise
        stylerclatex  = [];
        stylerclatex2 = [];
end

% precision od numbers
if nargin < 2 || isempty(precision)
    precision = 2;
end
if nargin > 1 && ~isempty(precision) && isscalar(precision)
    precision = round(precision);
end

% check on mandatory parameter
if nargin < 1
    error('Give me a table!');
end

% If T is an array, convert to a table
if isa(T,'numeric')
    T = array2table(T);
end

% if T is a timetable, convert to simple table
if istimetable(T)
    T = timetable2table(T);
end

% table dimensions
n_col     = size(T,2);
n_row     = size(T,1);

% If RowName is missing, add a ficticious one
if ~ismember('RowNames',T.Properties.VariableNames) && isempty(T.Properties.RowNames)
    T.Properties.RowNames = cellstr(string(1:size(T,1))');
    %T.Properties.RowNames = num2str((1:size(T,1))');
    %T.Properties.RowNames = matlab.lang.makeValidName(string((1:size(T,1))'));
end

% get row/column names
row_names = T.Properties.RowNames;
col_names = T.Properties.VariableNames;

% latex string initialization: add \bf and \it styles 
col_specs_latex  = repmat('r',1,n_col);
if ~isempty(stylerclatex)
    col_names_latex  = strjoin(col_names, [' & ' stylerclatex]); 
    col_names_latex  = [stylerclatex2 , col_names_latex];  % note: \bf vs \\bf
else
    col_names_latex  = strjoin(col_names, ' & ');
end

% disp string initialization
col_names_disp   = sprintf('%10s  | ', ' RAW NAME  ', col_names{:});

% latex string: arrange latex reserved characters/words
row_names_latex = strrep(row_names,'_','\_');
col_names_latex = strrep(col_names_latex,'_','\_');

% latex string: format of the table
if ~isempty(row_names_latex)
    col_specs_latex  = ['r'  col_specs_latex];
    col_names_latex = [' & ' col_names_latex];
end

% LaTeX string for file -- initialise file with header and tabular part.
fID = -1;
if ~isempty(filename) && ~strcmp(filename,'tabledisp_excel')
    while fID < 0
        [fID,~] = fopen(filename, 'w');
    end
    fprintf(fID, '\\begin{tabular}{%s}\n', col_specs_latex);
    fprintf(fID, '%s \\\\ \n', col_names_latex);
    fprintf(fID, '\\hline \n');
    fprintf(fID, '\\hline \n');
end

% LaTeX string for annotation/uitable -- \n are not used, as the
% interpreter wants to see a single line.
latex_string = '\begin{tabular}{';
latex_string = horzcat(latex_string , sprintf('%s} ', col_specs_latex));
latex_string = horzcat(latex_string , sprintf('%s \\\\ ', col_names_latex));
latex_string = horzcat(latex_string , sprintf('\\hline '));
latex_string = horzcat(latex_string , sprintf('\\hline '));

latex_string_full = latex_string;

% display string for the command window
%disp_string = horzcat('\n\n' , sprintf('%s \n\n ', col_names_disp));
disp_string = sprintf('%s \n\n ', col_names_disp);

% Add the table rows strings
try
    ns_col = n_col; start_col = 1;
    if ~isempty(row_names_latex)
        ns_col    = n_col + 1;
        start_col = 2;
    end
    for row = 1:n_row
        tempstring = cell(1,ns_col);
        tempstring(:) = {''};
        for col = 1:n_col
            value = T{row,col}; % content of a table cell
            while iscell(value), value = value{1,1}; end
            if ~isstring(value) && (isscalar(value) && isinf(value))
                value = '$\infty$';
            end
            if isstring(value) || ischar(value)
                tempstring{1,col+start_col-1} = char(value);
            elseif isdatetime(value)
                tempstring{1,col+start_col-1} = datestr(value); %#ok<DATST>
            else
                tempnumber = digitGouping(round(value,precision));
                tempstring{1,col+start_col-1} = num2str(tempnumber);
            end
        end
        tempstring_disp = tempstring;
        if ~isempty(row_names_latex)
            if ~isempty(stylerclatex)
                tempstring(1,1) = {strjoin({stylerclatex2 , row_names_latex{row}} , {' '})};
            else
                tempstring(1,1) = row_names_latex(row);
            end
            tempstring_disp(1,1) = row_names_latex(row);
        end

        tempstring_latex_sep      = strjoin(tempstring, ' & ');
        tempstring_latex_sep_full = strjoin(tempstring, ' & ');
        %tempstring_disp_sep = strjoin(tempstring, ' | ');
        tempstring_disp_sep  = sprintf('%10s  | ', tempstring_disp{:});
        %temparray = cellfun(@str2double,tempstring);
        %tempstring_disp_sep  = sprintf('% -10f  | ', temparray);

        % Long tables: estimate how many rows can say within 1200 characters
        if row==1
            lrow = numel(tempstring_latex_sep);
        else
            lrow = mean(lrow,numel(tempstring_latex_sep));
        end
        rmax = min(floor(1000/lrow),n_row);

        if rmax >= n_row
            % add the new line to the latex string
            latex_string = horzcat(latex_string , sprintf('%s \\\\ ', tempstring_latex_sep)); %#ok<AGROW>
        else
            if row < rmax - 4 || row == n_row - 1 || row == n_row
                % add the new line to the latex string
                latex_string = horzcat(latex_string , sprintf('%s \\\\ ', tempstring_latex_sep)); %#ok<AGROW>
            end                
            if (row == rmax - 4 || row == rmax - 3) 
                dummy = strjoin(repmat({' : '} , 1, n_col+1), ' & ') ; % repmat('& x ' , 1, n_col);
                latex_string = horzcat(latex_string , sprintf('%s \\\\ ', dummy)); %#ok<AGROW>
            end
        end

        % add the new line to the display string
        disp_string = horzcat(disp_string , sprintf('%s \n ', tempstring_disp_sep)); %#ok<AGROW>

        % create the full string in case of long tables, to be used in LaTeX documents
        latex_string_full = horzcat(latex_string_full , sprintf('%s \\\\ ', tempstring_latex_sep_full)); %#ok<AGROW> 

        if fID >= 0
            % add the new line to the file
            fprintf(fID, '%s \\\\ \n', tempstring_latex_sep_full);
        end

        clear temp;

    end
catch
    error('Something went wrong. Note that the table should only contain chars, strings or scalars.');
end

% Close the LaTeX string
latex_string = horzcat(latex_string , sprintf('\\hline '));
latex_string = horzcat(latex_string , sprintf('\\end{tabular}'));

% Close the LaTeX string_full
latex_string_full = horzcat(latex_string_full , sprintf('\\hline '));
latex_string_full = horzcat(latex_string_full , sprintf('\\end{tabular}'));

% To display the string in the command window, type:
% fprintf(disp_string);

% Close the LaTeX file
if fID >= 0
    fprintf(fID, '\\hline \n');
    fprintf(fID, '\\end{tabular}');
    fclose(fID);
end

% Export table also to excel
if strcmp(filename , 'tabledisp_excel')
    mfilename = [filename '.xlsx'];
    writetable(T,mfilename, 'WriteRowNames',true);
end


%% Apply comma as thousands separator.
    function formatted_number = digitGouping(number)
        [intPart , decPart]   = strtok(num2str(number),'.');
        ordermag = floor(log(abs(number))./log(10));
        intPart = intPart(end:-1:1);
        intPart = [sscanf(intPart,'%c',[3,inf])' repmat(',',ceil(length(intPart)/3),1)]';
        intPart = intPart(:)';
        intPart = deblank(intPart(1:(end-1)));
        if ordermag == 2
            formatted_number = strrep([intPart(end:-1:1) decPart],',','');
        else
            formatted_number = [intPart(end:-1:1) decPart];
        end
    end

end
%FScategory:UTIGEN


