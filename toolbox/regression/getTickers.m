function T = getTickers(varargin)
%getTickers retrieves representative or full ticker lists for a market/index
%
%<a href="matlab: docsearchFS('getTickers')">Link to the help function</a>
%
% Required input arguments:
%
%
% Optional input arguments:
%
%
% market :      Market or index name. Character or string scalar.
%               Possible values are:
%               'Nasdaq'   or 'Nasdaq100'
%               'SP500'    or 'S&P500'
%               'NYSE'
%               'Milan'    or 'FTSEMIB'
%               'London'   or 'FTSE100'
%               Default is 'Nasdaq'.
%               Example - 'market', 'Milan'
%               Data Types - char | string
%
% nStocks :     Number of representative stocks to return, excluding the
%               market index. Positive integer scalar.
%               Default is 10.
%               Maximum allowed value is 30 when Source is 'static'. When
%               Source is 'dynamic' Inf means all components.
%               Example - 'nStocks',15
%               Data Types - double
%
% Source :      Source used to build the market universe. Char or string.
%               Possible values are:
%               'static'  : use internal curated lists made up to 30
%               stocks.
%               'dynamic' : use dynamic download by scraping an online table of
%               index constituents, when available.
%               Default is 'static'.
%               Example - 'Source','dynamic'
%               Data Types - char | string
%
% RankByCap :   Use ranking criterion based on marketCap. Scalar boolean.
%               Retrieve market capitalization for each stock and rank the
%               constituents accordingly. Note that Ranking is applied only
%               to the constituents, not to the index row.
%               The default value of RankByCap is false.
%               Example - 'RankByCap',false
%               Data Types - logical
%
% IncludeIndex : Include the market index in the first row. Logical scalar.
%               Default is true.
%               If true, the first row of the output table contains the
%               market index and its descriptive name.
%               Example - 'IncludeIndex',false
%               Data Types - logical
%
% msg :         Display progress messages. Logical scalar.
%               Default is true.
%               Example - 'msg',false
%               Data Types - logical
%
% Output:
%
% T :           Output table.
%               Table with variables:
%               First column   = ticker symbol (ticker);
%               Second column  = company name (Name);
%               Third column = market capitalization (MarketCap) 
%                   This variable is present only if RankByCap=true.
%               The first row contains the market index if IncludeIndex=true.
%
% See also: getYahoo, getFundamentals
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('getTickers')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Default market and default number of stocks.
    % Get a table with two columns. Tickers and name from Nasdaq.
    % Use static representative list.
    T = getTickers;
    disp(T)
%}

%{
    % Top 10 constituents for Milan.
    % Use static representative list.
    T = getTickers('market','Milan');
    disp('Tickers from Milan stock exchange market')
    disp(T)
%}


%{
    % Top 10 constituents for S&P 500  without the index row.
    % Use static representative list.
    T = getTickers('market','SP500','IncludeIndex',false);
%}


%{
    % Top 15 Nasdaq tickers ranked by market cap using static list.
    T = getTickers('market','Nasdaq','nStocks',15,'Source','static');
%}

%{
    % Dynamic list then keep only top 20 by market cap.
    T = getTickers('market','London','nStocks',20,'Source','dynamic','RankByCap',true);
%}

%% Beginning of code


market = 'Nasdaq';
Source = 'static';
RankByCap = false;
IncludeIndex = true;
msg = true;
nStocks=10;

if nargin > 0
    options = struct('market',market,'Source',Source,'RankByCap',RankByCap, ...
        'IncludeIndex',IncludeIndex,'msg',msg,'nStocks',nStocks);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));

    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:getTickers:WrongInputOpt', ...
                'Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        aux.chkoptions(options,UserOptions);
    end

    for i = 1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end

    market = char(string(options.market));
    nStocks=options.nStocks;
    Source = char(string(options.Source));
    RankByCap=options.RankByCap;
    IncludeIndex = options.IncludeIndex;
    msg = options.msg;

end


if ~isscalar(nStocks) || ~isnumeric(nStocks) || isnan(nStocks) || nStocks < 1
    error('FSDA:getTickers:WrongInputOpt', ...
        'nStocks must be a positive integer scalar or Inf.');
end

if isfinite(nStocks) && floor(nStocks) ~= nStocks
    error('FSDA:getTickers:WrongInputOpt', ...
        'nStocks must be a positive integer scalar or Inf.');
end

if strcmp(Source,'static') && isfinite(nStocks) && nStocks > 30
    warning('FSDA:getTickers:WrongInputOpt', ...
        'Maximum allowed number of components is 30 when Source is ''static''. Using 30.');
    nStocks = 30;
end

Source = lower(strtrim(Source));
marketLower = lower(strtrim(market));

if isempty(Source)
    Source = 'static';
end

validSources = {'static','dynamic'};
if ~ismember(Source, validSources)
    error('FSDA:getTickers:WrongInputOpt', ...
        'Source must be ''static'' or ''dynamic''.');
end

%% Retrieve market universe (static or dynamic way)
% The output of localMarketUniverse is table U with two variables:
% ticker and Name.
% If RankByCap=true, variable MarketCap is added later.
% Depending on input option Source this information
% Can be retrieved from a curated list or in a dynamic way.
[indexSymbol, indexName, U] = localMarketUniverse(marketLower, Source);

if isempty(U) || height(U)==0
    error('FSDA:getTickers:NoUniverse', ...
        'Unable to retrieve constituents for market %s.', market);
end

%% Retrieving market capitalizations for the constituents
if RankByCap == true

    % Add ranking variable
    if msg == true
        fprintf('Retrieving market capitalizations for %d constituents of %s ...\n', ...
            height(U), market);
    end

    % Rank constituents
    U.MarketCap = localGetMarketCaps(U.ticker, msg);

    [~,ord] = sort(U.MarketCap, 'descend', 'MissingPlacement','last');
    U = U(ord,:);
end

%% After ordering take the first nStocks components
if isfinite(nStocks) && nStocks < height(U)
    U=U(1:nStocks,:);
end

%% Add market index row if requested
if IncludeIndex  == true
    if RankByCap == true
        Tidx = table(string(indexSymbol), string(indexName), NaN, ...
            'VariableNames', {'ticker','Name','MarketCap'});
    else
        Tidx = table(string(indexSymbol), string(indexName), ...
            'VariableNames', {'ticker','Name'});
    end
    T = [Tidx; U];
else
    T = U;
end

end

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [indexSymbol, indexName, U] = localMarketUniverse(marketLower, Source)

switch marketLower
    case {'nasdaq','nasdaq us','us tech','nasdaq100','nasdaq-100'}
        indexSymbol = "^IXIC";
        indexName   = "Nasdaq Composite";
        if strcmp(Source,'dynamic')
            url = 'https://en.wikipedia.org/wiki/Nasdaq-100';
            U = localWikiConstituents(url);
        else
            U = localStaticNasdaq();
        end

    case {'sp500','s&p500','sp-500','spx'}
        indexSymbol = "^GSPC";
        indexName   = "S&P 500 Index";
        if strcmp(Source,'dynamic')
            url = 'https://en.wikipedia.org/wiki/List_of_S%26P_500_companies';
            U = localWikiConstituents(url);
        else
            U = localStaticSP500();
        end

    case {'nyse','new york stock exchange'}
        indexSymbol = "^NYA";
        indexName   = "NYSE Composite";
        % For NYSE use static list in both cases unless you later add a
        % reliable dynamic source.
        U = localStaticNYSE();

    case {'milan','borsa italiana','italy','ftsemib'}
        indexSymbol = "FTSEMIB.MI";
        indexName   = "FTSE MIB Index";
        if strcmp(Source,'dynamic')
            url = 'https://en.wikipedia.org/wiki/FTSE_MIB';
            U = localWikiConstituents(url);
        else
            U = localStaticMilan();
        end

    case {'london','lse','london stock exchange','uk','ftse100'}
        indexSymbol = "^FTSE";
        indexName   = "FTSE 100 Index";
        if strcmp(Source,'dynamic')
            url = 'https://en.wikipedia.org/wiki/FTSE_100_Index';
            U = localWikiConstituents(url);
        else
            U = localStaticLondon();
        end

    otherwise
        error('FSDA:getTickers:WrongInputOpt', ...
            'Unknown market. Allowed values are ''Nasdaq'', ''SP500'', ''NYSE'', ''Milan'', ''London''.');
end

end

% -------------------------------------------------------------------------
function U = localWikiConstituents(url)

html = webread(url);
tree = htmlTree(html);
tables = findElement(tree,"table");

target = [];

for i = 1:numel(tables)
    txt = extractHTMLText(tables(i));

    if (contains(txt,"Ticker") || contains(txt,"Symbol") || contains(txt,"EPIC")) && ...
            (contains(txt,"Company") || contains(txt,"Security") || contains(txt,"Name"))
        target = tables(i);
        break
    end
end

if isempty(target)
    error('FSDA:getTickers:TableNotFound', ...
        'Unable to locate components table.');
end

rows = findElement(target,"tr");
if numel(rows) < 2
    error('FSDA:getTickers:BadTable', ...
        'The selected table does not contain enough rows.');
end

% Read header row dynamically
hdrCells = findElement(rows(1),"th");
if isempty(hdrCells)
    hdrCells = findElement(rows(1),"td");
end

if isempty(hdrCells)
    error('FSDA:getTickers:BadTable', ...
        'Unable to identify header cells in the constituents table.');
end

nHdr = numel(hdrCells);
hdrTxt = strings(nHdr,1);

for j = 1:nHdr
    hdrTxt(j) = lower(strtrim(extractHTMLText(hdrCells(j))));
end

idxTicker = find(contains(hdrTxt,"ticker") | contains(hdrTxt,"symbol") | contains(hdrTxt,"epic"), 1, 'first');
idxName   = find(contains(hdrTxt,"company") | contains(hdrTxt,"security") | contains(hdrTxt,"name"), 1, 'first');

if isempty(idxTicker)
    error('FSDA:getTickers:BadTable', ...
        'Unable to identify ticker column.');
end

if isempty(idxName)
    error('FSDA:getTickers:BadTable', ...
        'Unable to identify company/name column.');
end

ticker = strings(numel(rows)-1,1);
Name   = strings(numel(rows)-1,1);

k = 0;
for r = 2:numel(rows)
    cells = findElement(rows(r),"td");

    if numel(cells) >= max(idxTicker,idxName)
        sym = strtrim(extractHTMLText(cells(idxTicker)));
        nam = strtrim(extractHTMLText(cells(idxName)));

        if sym ~= "" && nam ~= ""
            k = k + 1;
            ticker(k) = localNormalizeTicker(sym,url);
            Name(k)   = nam;
        end
    end
end

ticker = ticker(1:k);
Name   = Name(1:k);

U = table(ticker, Name);

end

% -------------------------------------------------------------------------
function sym = localNormalizeTicker(sym,url)
sym = string(sym);
sym = strip(sym);
sym = replace(sym, ".", "-");

if contains(url,'FTSE_100_Index')
    if ~endsWith(sym,'.L')
        sym = sym + ".L";
    end
elseif contains(url,'Nasdaq-100')
    % keep US tickers as they are
elseif contains(url,'List_of_S%26P_500_companies')
    % keep US tickers as they are
end
end

% -------------------------------------------------------------------------
function U = localStaticNasdaq()

ticker = [ ...
    "AAPL"
    "MSFT"
    "NVDA"
    "AMZN"
    "GOOGL"
    "META"
    "TSLA"
    "AVGO"
    "COST"
    "NFLX"
    "ADBE"
    "AMD"
    "INTC"
    "CSCO"
    "PEP"
    "QCOM"
    "TMUS"
    "AMGN"
    "INTU"
    "TXN"
    "PYPL"
    "AMAT"
    "BKNG"
    "ISRG"
    "ADI"
    "MU"
    "LRCX"
    "PANW"
    "VRTX"
    "MELI"
    ];

Name = [ ...
    "Apple"
    "Microsoft"
    "NVIDIA"
    "Amazon"
    "Alphabet Class A"
    "Meta Platforms"
    "Tesla"
    "Broadcom"
    "Costco"
    "Netflix"
    "Adobe"
    "AMD"
    "Intel"
    "Cisco"
    "PepsiCo"
    "Qualcomm"
    "T-Mobile US"
    "Amgen"
    "Intuit"
    "Texas Instruments"
    "PayPal"
    "Applied Materials"
    "Booking Holdings"
    "Intuitive Surgical"
    "Analog Devices"
    "Micron Technology"
    "Lam Research"
    "Palo Alto Networks"
    "Vertex Pharmaceuticals"
    "MercadoLibre"
    ];

U = table(ticker, Name);

end

% -------------------------------------------------------------------------
function U = localStaticSP500()

ticker = [ ...
    "AAPL"
    "MSFT"
    "NVDA"
    "AMZN"
    "GOOGL"
    "META"
    "BRK-B"
    "TSLA"
    "UNH"
    "JNJ"
    "V"
    "XOM"
    "PG"
    "MA"
    "HD"
    "CVX"
    "MRK"
    "ABBV"
    "KO"
    "PEP"
    "COST"
    "AVGO"
    "ADBE"
    "CSCO"
    "WMT"
    "DIS"
    "MCD"
    "NFLX"
    "CRM"
    "BAC"
    ];

Name = [ ...
    "Apple"
    "Microsoft"
    "NVIDIA"
    "Amazon"
    "Alphabet Class A"
    "Meta Platforms"
    "Berkshire Hathaway B"
    "Tesla"
    "UnitedHealth"
    "Johnson & Johnson"
    "Visa"
    "Exxon Mobil"
    "Procter & Gamble"
    "Mastercard"
    "Home Depot"
    "Chevron"
    "Merck"
    "AbbVie"
    "Coca-Cola"
    "PepsiCo"
    "Costco"
    "Broadcom"
    "Adobe"
    "Cisco"
    "Walmart"
    "Disney"
    "McDonald's"
    "Netflix"
    "Salesforce"
    "Bank of America"
    ];

U = table(ticker, Name);

end

% -------------------------------------------------------------------------
function U = localStaticNYSE()

ticker = [ ...
    "BRK-B"
    "JPM"
    "V"
    "WMT"
    "XOM"
    "MA"
    "JNJ"
    "PG"
    "HD"
    "CVX"
    "KO"
    "BAC"
    "ABBV"
    "MRK"
    "MCD"
    "DIS"
    "IBM"
    "GE"
    "CAT"
    "GS"
    "MS"
    "AXP"
    "UNH"
    "T"
    "VZ"
    "RTX"
    "HON"
    "LOW"
    "BLK"
    "SCHW"
    ];

Name = [ ...
    "Berkshire Hathaway B"
    "JPMorgan Chase"
    "Visa"
    "Walmart"
    "Exxon Mobil"
    "Mastercard"
    "Johnson & Johnson"
    "Procter & Gamble"
    "Home Depot"
    "Chevron"
    "Coca-Cola"
    "Bank of America"
    "AbbVie"
    "Merck"
    "McDonald's"
    "Walt Disney"
    "IBM"
    "GE Aerospace"
    "Caterpillar"
    "Goldman Sachs"
    "Morgan Stanley"
    "American Express"
    "UnitedHealth"
    "AT&T"
    "Verizon"
    "RTX"
    "Honeywell"
    "Lowe's"
    "BlackRock"
    "Charles Schwab"
    ];

U = table(ticker, Name);

end

% -------------------------------------------------------------------------
function U = localStaticMilan()
% NOTE:
% This is a curated list of major Italian stocks (approx. FTSE MIB).
% It is not guaranteed to exactly match the current index composition.
ticker = [ ...
    "ENEL.MI"
    "ENI.MI"
    "ISP.MI"
    "UCG.MI"
    "G.MI"
    "RACE.MI"
    "STMMI.MI"
    "TIT.MI"
    "SRG.MI"
    "PST.MI"
    "MB.MI"
    "MONC.MI"
    "LDO.MI"
    "CPR.MI"
    "BZU.MI"
    "BAMI.MI"
    "TEN.MI"
    "HER.MI"
    "A2A.MI"
    "INW.MI"
    "NEXI.MI"
    "ERG.MI"
    "ACE.MI"
    "DIA.MI"
    "AMP.MI"
    "FBK.MI"
    "TRN.MI"
    "PRY.MI"
    "BPE.MI"
    "IVG.MI"
    ];

Name = [ ...
    "Enel"
    "Eni"
    "Intesa Sanpaolo"
    "UniCredit"
    "Assicurazioni Generali"
    "Ferrari"
    "STMicroelectronics"
    "Telecom Italia"
    "Snam"
    "Poste Italiane"
    "Mediobanca"
    "Moncler"
    "Leonardo"
    "Campari"
    "Buzzi"
    "Banco BPM"
    "Tenaris"
    "Hera"
    "A2A"
    "Inwit"
    "Nexi"
    "ERG"
    "Acea"
    "DiaSorin"
    "Amplifon"
    "FinecoBank"
    "Terna"
    "Prysmian"
    "BPER Banca"
    "Iveco Group"
    ];

U = table(ticker, Name);

end

% -------------------------------------------------------------------------
function U = localStaticLondon()

ticker = [ ...
    "SHEL.L"
    "AZN.L"
    "HSBA.L"
    "ULVR.L"
    "BP.L"
    "GSK.L"
    "RIO.L"
    "DGE.L"
    "NG.L"
    "BARC.L"
    "LSEG.L"
    "LLOY.L"
    "REL.L"
    "GLEN.L"
    "VOD.L"
    "PRU.L"
    "AAL.L"
    "SSE.L"
    "BA.L"
    "BT-A.L"
    "IMB.L"
    "STAN.L"
    "CRH.L"
    "SMIN.L"
    "NWG.L"
    "RR.L"
    "CNA.L"
    "ABF.L"
    "AHT.L"
    "ADM.L"
    ];

Name = [ ...
    "Shell"
    "AstraZeneca"
    "HSBC"
    "Unilever"
    "BP"
    "GSK"
    "Rio Tinto"
    "Diageo"
    "National Grid"
    "Barclays"
    "London Stock Exchange Group"
    "Lloyds Banking Group"
    "RELX"
    "Glencore"
    "Vodafone"
    "Prudential"
    "Anglo American"
    "SSE"
    "BAE Systems"
    "BT Group"
    "Imperial Brands"
    "Standard Chartered"
    "CRH"
    "Smiths Group"
    "NatWest"
    "Rolls-Royce"
    "Centrica"
    "Associated British Foods"
    "Ashtead Group"
    "Admiral Group"
    ];

U = table(ticker, Name);

end

% -------------------------------------------------------------------------
function mc = localGetMarketCaps(symbols, msg)

n = numel(symbols);
mc = nan(n,1);

try
    F = getFundamentals(symbols,'Fields',{'TickerSymbol','marketCap'});
    if ismember('TickerSymbol', F.Properties.VariableNames) && ...
            ismember('marketCap', F.Properties.VariableNames)

        fsym = string(F.TickerSymbol);
        fmc  = F.marketCap;

        for i = 1:n
            idx = find(fsym == string(symbols(i)),1,'first');
            if ~isempty(idx)
                if isnumeric(fmc)
                    mc(i) = fmc(idx);
                elseif iscell(fmc)
                    if isempty(fmc{idx})
                        mc(i) = NaN;
                    elseif isnumeric(fmc{idx}) || islogical(fmc{idx})
                        mc(i) = double(fmc{idx});
                    end
                end
            end
        end
        return
    end
catch
    if msg
        fprintf('getFundamentals unavailable or failed. Falling back to individual retrieval.\n');
    end
end

for i = 1:n
    try
        mc(i) = getMarketCap(char(symbols(i)));
    catch
        mc(i) = NaN;
    end
end

end

