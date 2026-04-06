function T = getTickers(varargin)
%getTickers retrieves representative or full ticker lists for a market/index
%
%<a href="matlab: docsearchFS('getTickers')">Link to the help function</a>
%
%   getTickers, getYahoo and getFundamentals can be used jointly to build a
%   complete workflow: from the selection of representative market tickers,
%   to the retrieval and dynamic interactive plot of their price time
%   series, and finally to the extraction of their fundamental financial
%   information.
%   For background on financial data and market analysis, see:
%   Yahoo Finance API documentation https://finance.yahoo.com/ 
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
%               'DAX'      or 'DAX40'
%               'CAC40'    or 'CAC'
%               'Nikkei225' or 'Nikkei'
%               Default is 'Nasdaq'.
%               Example - 'market', 'DAX'
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
% See also: getYahoo, getFundamentals, rsindex, candle, movavg
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

%{
    % Top 10 DAX constituents using static representative list.
    T = getTickers('market','DAX');
%}

%{
    % Dynamic CAC 40 constituents.
    T = getTickers('market','CAC40','Source','dynamic','nStocks',Inf);
%}

%{
    % Top 20 Nikkei names ranked by market cap.
    T = getTickers('market','Nikkei225','nStocks',20,'RankByCap',true);
%}

%{
    Example of combined use of getTickers and getFundamentals.
    T = getTickers('market','SP500','nStocks',10,'RankByCap',true);
    disp(T)

    % Retrieve fundamentals and verify ranking
    F = getFundamentals(T.ticker(2:end),'Fields','basic');
    disp(F(:,{'TickerSymbol','marketCap'}))
%}

%{
    % Example of combined use of getTickers, getYahoo and getFundamentals.
    % Step 1: get tickers
    T = getTickers('market','CAC40','nStocks',5,'RankByCap',true);

    % Step 2: get prices
    P = getYahoo(T.ticker(2:end));

    % Step 3: fundamentals
    F = getFundamentals(T.ticker(2:end),'Fields','performance');

    % Combine outputs manually
    disp(T)
    disp(F)
%}

%{
    %% Example with full pipeline.
    % STEP 1: Retrieve top tickers (static, no ranking).
    T = getTickers('market','Nasdaq','nStocks',5);
    disp(T)

    % STEP 2: Download price data for those tickers.
    out = getYahoo(T.ticker(2:end));  % skip index
    disp(out)

    % STEP 3: Retrieve fundamentals.
    F = getFundamentals(T.ticker(2:end),'Fields','basic');
    disp(F)
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

if isempty(Source)
    Source = 'static';
end

Source = lower(strtrim(Source));

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

marketLower = lower(strtrim(market));


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
% can be retrieved from a curated list or in a dynamic way.
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

    case {'dax','germany','de','dax40'}
        indexSymbol = "^GDAXI";
        indexName   = "DAX Index";
        if strcmp(Source,'dynamic')
            url = 'https://en.wikipedia.org/wiki/DAX';
            U = localWikiConstituents(url);
        else
            U = localStaticDAX();
        end

    case {'cac','cac40','france','paris'}
        indexSymbol = "^FCHI";
        indexName   = "CAC 40 Index";
        if strcmp(Source,'dynamic')
            url = 'https://en.wikipedia.org/wiki/CAC_40';
            U = localWikiConstituents(url);
        else
            U = localStaticCAC40();
        end

    case {'nikkei','nikkei225','japan','jp'}
        indexSymbol = "^N225";
        indexName   = "Nikkei 225 Index";
        % Dynamic Wikipedia parsing is not straightforward here because the
        % page uses sector sublists rather than a simple ticker table.

        if strcmp(Source,'dynamic')
            warning('FSDA:getTickers:NoDynamicSource', ...
                ['Dynamic retrieval for Nikkei 225 is currently unavailable. ' ...
                'Using static curated list instead.']);
        end
        U = localStaticNikkei225();

    otherwise
        error('FSDA:getTickers:WrongInputOpt', ...
            ['Unknown market. Allowed values are ''Nasdaq'', ''SP500'', ''NYSE'', ' ...
            '''Milan'', ''London'', ''DAX'', ''CAC40'', ''Nikkei225''.']);
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

elseif contains(url,'DAX')
    if ~endsWith(sym,'.DE')
        sym = sym + ".DE";
    end

elseif contains(url,'CAC_40')
    if ~endsWith(sym,'.PA') && ~endsWith(sym,'.AS')
        % Most CAC members are .PA, but ArcelorMittal is commonly MT.AS
        if sym == "MT"
            sym = "MT.AS";
        else
            sym = sym + ".PA";
        end
    end

    % elseif contains(url,'Nikkei_225')
    %     if ~endsWith(sym,'.T')
    %         sym = sym + ".T";
    %     end

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
function U = localStaticDAX()

ticker = [ ...
    "SAP.DE"
    "SIE.DE"
    "ALV.DE"
    "AIR.DE"
    "DTE.DE"
    "MUV2.DE"
    "MBG.DE"
    "RHM.DE"
    "BMW.DE"
    "BAS.DE"
    "BAYN.DE"
    "ADS.DE"
    "IFX.DE"
    "DB1.DE"
    "VOW3.DE"
    "HEN3.DE"
    "SHL.DE"
    "DHL.DE"
    "EOAN.DE"
    "ENR.DE"
    "BEI.DE"
    "HEI.DE"
    "QIA.DE"
    "FRE.DE"
    "CON.DE"
    "HNR1.DE"
    "MRK.DE"
    "MTX.DE"
    "SY1.DE"
    "P911.DE"
    ];

Name = [ ...
    "SAP"
    "Siemens"
    "Allianz"
    "Airbus"
    "Deutsche Telekom"
    "Munich Re"
    "Mercedes-Benz Group"
    "Rheinmetall"
    "BMW"
    "BASF"
    "Bayer"
    "Adidas"
    "Infineon"
    "Deutsche Boerse"
    "Volkswagen Pref"
    "Henkel Pref"
    "Siemens Healthineers"
    "DHL Group"
    "E.ON"
    "Siemens Energy"
    "Beiersdorf"
    "Heidelberg Materials"
    "QIAGEN"
    "Fresenius"
    "Continental"
    "Hannover Rueck"
    "Merck KGaA"
    "MTU Aero Engines"
    "Symrise"
    "Porsche AG"
    ];

U = table(ticker, Name);

end

function U = localStaticCAC40()

ticker = [ ...
    "MC.PA"
    "OR.PA"
    "TTE.PA"
    "SAN.PA"
    "AIR.PA"
    "AI.PA"
    "SU.PA"
    "EL.PA"
    "BNP.PA"
    "CS.PA"
    "SAF.PA"
    "DG.PA"
    "RI.PA"
    "CAP.PA"
    "ORA.PA"
    "ACA.PA"
    "ENGI.PA"
    "VIE.PA"
    "WLN.PA"
    "HO.PA"
    "LR.PA"
    "KER.PA"
    "ML.PA"
    "CA.PA"
    "SGO.PA"
    "STM.PA"
    "BN.PA"
    "PUB.PA"
    "TEP.PA"
    "AC.PA"
    ];

Name = [ ...
    "LVMH"
    "L'Oreal"
    "TotalEnergies"
    "Sanofi"
    "Airbus"
    "Air Liquide"
    "Schneider Electric"
    "EssilorLuxottica"
    "BNP Paribas"
    "AXA"
    "Safran"
    "Vinci"
    "Pernod Ricard"
    "Capgemini"
    "Orange"
    "Credit Agricole"
    "Engie"
    "Veolia"
    "Worldline"
    "Thales"
    "Legrand"
    "Kering"
    "Michelin"
    "Carrefour"
    "Saint-Gobain"
    "STMicroelectronics"
    "Danone"
    "Publicis"
    "Teleperformance"
    "Accor"
    ];

U = table(ticker, Name);

end

function U = localStaticNikkei225()

ticker = [ ...
    "7203.T"
    "6758.T"
    "6861.T"
    "9983.T"
    "9984.T"
    "8306.T"
    "8035.T"
    "9432.T"
    "9433.T"
    "4063.T"
    "6501.T"
    "7974.T"
    "7267.T"
    "6954.T"
    "4519.T"
    "6098.T"
    "8766.T"
    "8058.T"
    "8031.T"
    "6981.T"
    "6367.T"
    "7751.T"
    "6857.T"
    "2914.T"
    "4502.T"
    "4568.T"
    "6971.T"
    "4661.T"
    "6902.T"
    "6503.T"
    ];

Name = [ ...
    "Toyota Motor"
    "Sony Group"
    "Keyence"
    "Fast Retailing"
    "SoftBank Group"
    "Mitsubishi UFJ Financial Group"
    "Tokyo Electron"
    "NTT"
    "KDDI"
    "Shin-Etsu Chemical"
    "Hitachi"
    "Nintendo"
    "Honda Motor"
    "Fanuc"
    "Chugai Pharmaceutical"
    "Recruit Holdings"
    "Tokio Marine Holdings"
    "Mitsubishi Corp"
    "Mitsui & Co"
    "Murata Manufacturing"
    "Daikin Industries"
    "Canon"
    "Advantest"
    "Japan Tobacco"
    "Takeda Pharmaceutical"
    "Daiichi Sankyo"
    "Kyocera"
    "Oriental Land"
    "Denso"
    "Mitsubishi Electric"
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

%FScategory:UTI-FIN