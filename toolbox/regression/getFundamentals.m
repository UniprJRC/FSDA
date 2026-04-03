function T = getFundamentals(ticker)
%getFundamentals Retrieve Yahoo Finance information for one or more tickers
%
%<a href="matlab: docsearchFS('getFundamentals')">Link to the help function</a>
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
% Output:
%
% T       :     Table with one row for each ticker.
%               The Rownames of this table contains the company name.
%               The first column contains the ticker symbol.
%               The remaining columns contain all available information
%               returned by Yahoo Finance through Python yfinance module.
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
        'symbols must be a char, string, string array, or cell array of char.');
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

    % Extract company name with fallback rules
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
    allFields = union(allFields, fieldnames(info));
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

% Convert to table
T = struct2table(S, 'AsArray', true);

% Move FirmName and TickerSymbol to the first two columns
vnames = T.Properties.VariableNames;

firstVars = {'CompanyName','TickerSymbol'};
remaining = setdiff(vnames, firstVars, 'stable');
T = T(:, [firstVars remaining]);
T.Properties.RowNames=T{:,1};
T=T(:,2:end);
end

function infoStruct = getFromPythonModule_yfinance(symbol)
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
    ticker = yf.Ticker(symbol);
    infoPy = ticker.info;

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
                out = char(val);
            catch
                out = val;
            end
        end
    end
end