function T = getFundamentals(ticker, varargin)
%getFundamentals retrieves detailed company metrics for supplied tickers. 
%
%<a href="matlab: docsearchFS('getFundamentals')">Link to the help function</a>
%
%   getTickers, getYahoo and getFundamentals can be used jointly to build a
%   complete workflow: from the selection of representative market tickers,
%   to the retrieval and dynamic interactive plot of their price time
%   series, and finally to the extraction of their fundamental financial
%   information.
%   For background on financial data and market analysis, see:
%   Yahoo Finance API documentation https://finance.yahoo.com/ 
%   
%
% Required input arguments:
%
% ticker :     Ticker symbol(s). Character, string, string array or cell array of char.
%               It can be:
%               - a character vector, for example 'ENEL.MI';
%               - a string scalar, for example "ENEL.MI";
%               - a string array, for example ["ENEL.MI" "AAPL" "MSFT"];
%               - a cell array of char, for example {'ENEL.MI','AAPL','MSFT'};
%
% Optional input arguments:
%
% RankBy :     Rank rows by a numeric variable in descending order.
%              If the requested field is missing, an error is produced.
%              Missing numeric values are placed at the bottom.
%              Default is '', which means no ranking. Note that ranking
%              happens before field reduction
%              Example - 'RankBy','marketCap'
%              Data Types - char | string
%
% ExplainFields : Display explanation of available fields. Logical or string.
%               If true, prints a description of the most relevant fields.
%               If equal to:
%                   'basic'       → explain basic preset
%                   'valuation'   → explain valuation preset
%                   'performance' → explain performance preset
%               Default is false.
%               Example - 'ExplainFields','valuation'
%               Data Types - logical | char | string
%
% Fields : Select subset of output variables.
%              Possible values:
%              'all'         = return all available fields (default)
%
%              'basic'       = {'CompanyName','TickerSymbol',...
%                               'marketCap','sector','industry'}
%
%              'valuation'   = {'CompanyName','TickerSymbol',...
%                               'marketCap','enterpriseValue',...
%                               'trailingPE','forwardPE','priceToBook'}
%
%              'performance' = {'CompanyName','TickerSymbol',...
%                               'returnOnEquity','profitMargins','revenueGrowth'}
%
%              Alternatively, a custom list of fields can be supplied:
%              For instance 'Fields',{'marketCap','trailingPE'}
%              Example - 'Fields','basic'
%              Data Types - char | string | cell array of characters
% Output:
%
% T       :     Table with one row for each ticker. Table with k columns.
%               The company name is stored in the row names (making sure
%               that there are no repeated elements and that the names are valid)
%               and is repeated as the first table variable. The second
%               variable of the output table is TickerSymbol. The remaining
%               columns contain all available information returned by Yahoo
%               Finance through Python yfinance module.
%
%
% See also: getYahoo, getTickers, rsindex, candle, movavg
%
% References:
%
%  Damodaran, A. (2012). "Investment Valuation: Tools and Techniques for
%  Determining the Value of Any Asset", 3rd Edition, Wiley.
%
%  Cochrane, J. H. (2023). "Asset Pricing, Revised Edition",
%  Princeton University Press, Princeton.
%
%  Koller, T., Goedhart, M., and Wessels, D. (2020).
%  "Valuation: Measuring and Managing the Value of Companies,
%  7th Edition", Wiley, Hoboken.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('getFundamentals')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples:
%
%{
    % One ticker.
    T = getFundamentals('ENEL.MI');
%}

%{
    % Multiple tickers Example 1.
    % The input is a cell array of character vectors.
    T = getFundamentals({'ENEL.MI','AAPL','MSFT'});
%}

%{
    % Multiple tickers Example 2.
    % The input is a string array.
    T = getFundamentals(["ENEL.MI","RACE.MI","AAPL"]);
%}

%{
    % Basic information only.
    T = getFundamentals('AAPL','Fields','basic');
%}

%{
    % Valuation metrics.
    T = getFundamentals(["AAPL","MSFT","NVDA"],'Fields','valuation');
%}

%{
    % Performance metrics.
    T = getFundamentals('ENEL.MI','Fields','performance');
%}

%{
    % Custom fields.
    T = getFundamentals('AAPL','Fields',{'marketCap','trailingPE'});
%}

%{
    %% Rank by market capitalization
    T = getFundamentals(["AAPL","MSFT","NVDA"], ...
        'Fields','basic','RankBy','marketCap');
    disp(T)
%}

%{
    % Rank by forward PE
    T = getFundamentals(["AAPL","MSFT","NVDA"], ...
        'Fields','valuation','RankBy','forwardPE');
%}

%{
    % Rank performance metrics.
    T = getFundamentals(["AAPL","MSFT","NVDA"], ...
        'Fields','performance','RankBy','returnOnEquity');
%}

%{
    %% Example of combined use of getTickers with getFundamentals.
    T = getTickers('market','DAX','nStocks',8,'RankByCap',true);

    % Retrieve selected fields only
    F = getFundamentals(T.ticker(2:end), ...
        'Fields',{'marketCap','trailingPE','returnOnEquity'});

    disp(F)
%}



%% Beginning of code

if nargin < 1 || isempty(ticker)
    error('FSDA:getFundamentals:MissingInput', ...
        'At least one ticker symbol must be supplied.');
end

% Normalize input to string column vector
if ischar(ticker)
    ticker = string({ticker});
elseif isstring(ticker)
    ticker = ticker(:);
elseif iscell(ticker) && all(cellfun(@ischar, ticker))
    ticker = string(ticker(:));
else
    error('FSDA:getFundamentals:WrongInput', ...
        'ticker must be a char, string, string array, or cell array of char.');
end

% Ensure Python is available
pe = pyenv;
assert(~isempty(pe.Version), ...
    'Python is not configured in MATLAB. Run pyenv to configure.');

% Import yfinance
try
    py.importlib.import_module('yfinance');
catch
    error('FSDA:getFundamentals:MissingPackage', ...
        'Python package "yfinance" is not installed. Run: python.exe -m pip install yfinance');
end

n = numel(ticker);
Fields='all';
RankBy='';

ExplainFields = false;


if nargin>1

    options = struct('Fields',Fields,'RankBy',RankBy, 'ExplainFields',ExplainFields);
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:getFundamentals:WrongInputOpt', ...
                'Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        aux.chkoptions(options,UserOptions);
    end


    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
    Fields=options.Fields;
    RankBy = options.RankBy;
    ExplainFields=options.ExplainFields;
end

if ~isequal(ExplainFields,false)
    explainFieldsLocal(ExplainFields);
end


% Use cell storage first, because structs may have different fields
infoCell = cell(n,1);
allFields = {};

for i = 1:n
    try
        info = getFromPythonModule_yfinance(ticker(i));
    catch ME
        warning('FSDA:getFundamentals:DownloadFailed', ...
            'Unable to retrieve info for ticker %s. Error: %s', ...
            char(ticker(i)), ME.message);
        info = struct();
    end

    % Add ticker symbol explicitly
    info.TickerSymbol = char(ticker(i));

    % Extract company name with fallback rules (that is companyName is
    % taken from 'longName' or 'shortName' or 'displayName')
    if isfield(info,'longName') && ~isempty(info.longName)
        companyName = info.longName;
    elseif isfield(info,'shortName') && ~isempty(info.shortName)
        companyName = info.shortName;
    elseif isfield(info,'displayName') && ~isempty(info.displayName)
        companyName = info.displayName;
    else
        companyName = char(ticker(i));
    end
    info.CompanyName = companyName;

    infoCell{i} = info;
    %allFields = union(allFields, fieldnames(info));

    fnNew = fieldnames(info);
    % Note that each ticker symbol can have different variables
    % In every iteration we add the fields which are not present yet
    allFields = [allFields; setdiff(fnNew, allFields, 'stable')];
end

% Rebuild homogeneous struct array
allFields = allFields(:);
S = repmat(cell2struct(cell(numel(allFields),1), allFields, 1), n, 1);

for i = 1:n
    info = infoCell{i};
    fn = fieldnames(info);
    for j = 1:numel(fn)
        S(i).(fn{j}) = info.(fn{j});
    end
end

% ============================================================
% FORCE NUMERIC CONSISTENCY
% ============================================================

fn = fieldnames(S);

for j = 1:numel(fn)
    f = fn{j};

    col = {S.(f)}';   % extract column as cell

    % Check if column is numeric-like
    isNumLike = cellfun(@(x) isnumeric(x) || islogical(x) || isempty(x), col);

    if all(isNumLike)
        tmp = nan(numel(col),1);

        for i = 1:numel(col)
            if isempty(col{i})
                % Put NaN for missing numeric values
                % In this way Column stays double, and does not become a cell
                tmp(i) = NaN;   
            else
                tmp(i) = double(col{i});
            end
        end

        % overwrite struct field with numeric array
        for i = 1:numel(col)
            S(i).(f) = tmp(i);
        end
    end
end

T = struct2table(S, 'AsArray', true);

% ============================================================
% OPTIONAL RANKING
% ============================================================
if ~isempty(RankBy)
    RankBy = char(string(RankBy));

    if ~ismember(RankBy, T.Properties.VariableNames)
        error('FSDA:getFundamentals:WrongInputOpt', ...
            'Requested RankBy field "%s" is not present in the output table.', RankBy);
    end

    rankVar = T.(RankBy);

    % Convert ranking variable to numeric if possible
    if isnumeric(rankVar) || islogical(rankVar)
        rankNum = double(rankVar);
    elseif iscell(rankVar)
        rankNum = nan(height(T),1);
        for i = 1:height(T)
            try
                if isnumeric(rankVar{i}) || islogical(rankVar{i})
                    rankNum(i) = double(rankVar{i});
                elseif ischar(rankVar{i}) || isstring(rankVar{i})
                    tmp = str2double(string(rankVar{i}));
                    if ~isnan(tmp)
                        rankNum(i) = tmp;
                    end
                end
            catch
            end
        end
    elseif isstring(rankVar) || ischar(rankVar)
        rankNum = str2double(string(rankVar));
    else
        error('FSDA:getFundamentals:WrongInputOpt', ...
            'Field "%s" cannot be used for ranking because it is not numeric.', RankBy);
    end

    if all(isnan(rankNum))
    error('FSDA:getFundamentals:WrongInputOpt', ...
        'Field "%s" cannot be used for ranking because no numeric values could be extracted.', RankBy);
    end

    % Missing values go to the bottom
    rankNum(isnan(rankNum)) = -Inf;

    [~, ord] = sort(rankNum, 'descend');
    T = T(ord,:);
end

% ============================================================
% FIELD SELECTION LOGIC
% ============================================================

if isempty(Fields)
    Fields='all';
end


% Define presets
switch lower(string(Fields))

    case "basic"
        reqFields = {'CompanyName','TickerSymbol', ...
            'marketCap','sector','industry'};

    case "valuation"
        reqFields = {'CompanyName','TickerSymbol', ...
            'marketCap','enterpriseValue','trailingPE', ...
            'forwardPE','priceToBook'};

    case "performance"
        reqFields = {'CompanyName','TickerSymbol', ...
            'returnOnEquity','profitMargins','revenueGrowth'};

    case "all"
        reqFields=T.Properties.VariableNames;

    otherwise
        % Custom user fields
        if ischar(Fields)
            reqFields = {Fields};
        elseif isstring(Fields)
            reqFields = cellstr(Fields(:));
        elseif iscell(Fields)
            reqFields = Fields(:)';
        else
            error('FSDA:getFundamentals:WrongInputOpt', ...
                'Fields must be ''all'', a preset name, or a cell array.');
        end
end



% Keep only fields that exist
existing = T.Properties.VariableNames;
reqFields = intersect(reqFields, existing, 'stable');
% Always keep first two columns if available
baseVars = {'CompanyName','TickerSymbol'};
baseVars = intersect(baseVars, existing, 'stable');

T = T(:, unique([baseVars reqFields], 'stable'));

rn = string(T.CompanyName);
rn(ismissing(rn)) = "UnknownCompany";
rn(rn=="") = "UnknownCompany";
rn = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(rn));
T.Properties.RowNames = cellstr(rn);

% T=T(:,2:end);
end

function infoStruct = getFromPythonModule_yfinance(ticker)
% getFromPythonModule_yfinance  Retrieve ALL available information for a ticker from Yahoo Finance
% using Python's yfinance library as a MATLAB struct.
%
%   infoStruct = getFromPythonModule_yfinance('ENEL.MI')

% Ensure Python is available
assert(~isempty(pyenv().Version), ...
    'Python is not configured in MATLAB. Run pyenv to configure.');

% Import yfinance
yf = py.importlib.import_module('yfinance');

% Fetch ticker info (Python dict)
tobj  = yf.Ticker(ticker);
infoPy = tobj .info;

% Convert Python dict to MATLAB struct
infoStruct = pyDict2Matlab(infoPy);
end

% Convert python dictionary to MATALB
function out = pyDict2Matlab(pyDict)
out = struct();
keys = cell(py.list(pyDict.keys()));

for i = 1:numel(keys)
    keyPy = keys{i};
    key = matlab.lang.makeValidName(char(keyPy));
    val = pyDict{keyPy};
    out.(key) = convertPy(val);
end
end

function out = convertPy(val)

if isa(val, 'py.NoneType')
    out = [];

elseif isa(val, 'py.float')
    out = double(val);

elseif isa(val, 'py.int')
    out = double(val);

elseif isa(val, 'py.bool')
    out = logical(val);

elseif isa(val, 'py.str')
    out = char(val);

elseif isa(val, 'py.list') || isa(val, 'py.tuple')
    out = cellfun(@convertPy, cell(val), 'UniformOutput', false);

elseif isa(val, 'py.dict')
    out = pyDict2Matlab(val);

else
    try
        out = double(val);
    catch
        try
            out = logical(val);
        catch
            try
                out = char(val);
            catch
                out = string(py.str(val));
            end
        end
    end
end
end

function explainFieldsLocal(mode)

fprintf('\n============================================\n');
fprintf('AVAILABLE FIELDS IN getFundamentals\n');
fprintf('============================================\n');

basicFields = {
    'CompanyName        : Company name'
    'TickerSymbol       : Ticker symbol'
    'marketCap          : Market capitalization'
    'sector             : Economic sector'
    'industry           : Industry classification'
};

valuationFields = {
    'marketCap          : Market capitalization'
    'enterpriseValue    : Firm value (incl. debt)'
    'trailingPE         : Price/Earnings (past)'
    'forwardPE          : Price/Earnings (forecast)'
    'priceToBook        : Price / Book value'
};

performanceFields = {
    'returnOnEquity     : ROE'
    'profitMargins      : Net profit margin'
    'revenueGrowth      : Revenue growth'
};

otherFields = {
    'currentPrice       : Latest price'
    'beta               : Market risk'
    'dividendYield      : Dividend yield'
    'totalRevenue       : Total revenue'
    'grossMargins       : Gross margin'
    'operatingMargins   : Operating margin'
    'freeCashflow       : Free cash flow'
};

if islogical(mode)
    if mode
        mode = "all";
    else
        return
    end
else
    mode = lower(string(mode));
end

switch mode
    case "all"
        fprintf('\n--- BASIC ---\n');       fprintf('%s\n', basicFields{:});
        fprintf('\n--- VALUATION ---\n');   fprintf('%s\n', valuationFields{:});
        fprintf('\n--- PERFORMANCE ---\n'); fprintf('%s\n', performanceFields{:});
        fprintf('\n--- OTHER ---\n');       fprintf('%s\n', otherFields{:});

    case "basic"
        fprintf('\n--- BASIC ---\n'); fprintf('%s\n', basicFields{:});

    case "valuation"
        fprintf('\n--- VALUATION ---\n'); fprintf('%s\n', valuationFields{:});

    case "performance"
        fprintf('\n--- PERFORMANCE ---\n'); fprintf('%s\n', performanceFields{:});

    otherwise
        warning('FSDA:getFundamentals:WrongInputOpt', ...
            'Unknown ExplainFields option.');
end

fprintf('============================================\n\n');

end

%FScategory:UTI-FIN