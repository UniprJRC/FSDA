function T = getTickers(market, nStocks)
%getTickers returns representative ticker symbols for a market
%
%<a href="matlab: docsearchFS('getTickers')">Link to the help function</a>
%
%
% Required input arguments:
%
% Optional input arguments:
%
% market :      Market name. Character or string scalar.
%               Possible values are:
%               'Nasdaq' (default);
%               'NYSE';
%               'Milan';
%               'London'.
%               The first row of the output is always the overall market
%               index.
%               Example - 'Milan'
%               Data Types - char | string
%
% nStocks :     Number of representative stocks to return, excluding the
%               market index. Positive integer scalar.
%               Default is 10.
%               Maximum allowed value is 100.
%               Example - 15
%               Data Types - double
%
% Output:
%
% TT      :     Ticker and Name. Table.
%               Table with 2 variables. 
%               First column = ticker symbol
%               Second column company or index name.
%               The first row always contains the market index.
%
%
% See also: getYahoo
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
    T = getTickers;
    disp(T)
%}

%{
    % Ten representative tickers for Milan.
     T = getTickers('Milan');
     disp(T)
%}

%{
    % Fifteen representative tickers for Nasdaq.
     T = getTickers('Nasdaq',15);
%}

%{
    % Twenty representative tickers for NYSE.
     T = getTickers('NYSE',20);
%}
%
%{
% Use returned symbols with getYahoo.
    T = getTickers('London',5);
    out = getYahoo(T.Symbol(2:end));
%}

%% Beginning of code

if nargin < 1 || isempty(market)
    market = 'Nasdaq';
end

if nargin < 2 || isempty(nStocks)
    nStocks = 10;
end

market = char(string(market));

if ~isscalar(nStocks) || ~isnumeric(nStocks) || isnan(nStocks) || ...
        nStocks < 1 || floor(nStocks) ~= nStocks
    error('FSDA:getTickers:WrongInputOpt', ...
        'nStocks must be a positive integer scalar.');
end

if nStocks > 100
    error('FSDA:getTickers:WrongInputOpt', ...
        'nStocks cannot be greater than 100.');
end

marketLower = lower(strtrim(market));

switch marketLower
    case {'nasdaq','nasdaq us','us tech','nasdaq100','nasdaq-100'}
        Tfull = localNasdaq();

    case {'nyse','new york stock exchange'}
        Tfull = localNYSE();

    case {'milan','borsa italiana','italy','ftsemib'}
        Tfull = localMilan();

    case {'london','lse','london stock exchange','uk','ftse100'}
        Tfull = localLondon();

    case {'sp500','s&p500','sp-500','standard and poor','standard & poor'}
        Tfull = localSP500();

    otherwise
        error('FSDA:getTickers:WrongInputOpt', ...
            'Unknown market. Allowed markets are ''Nasdaq'', ''NYSE'', ''Milan'', ''London'', ''SP500''.');
end

% First row = market index, following rows = equities
nAvailableStocks = height(Tfull) - 1;
nToTake = min(nStocks, nAvailableStocks);

T = Tfull([1; (2:nToTake+1)'], :);

end

% -------------------------------------------------------------------------
function T = localNasdaq()

Symbol = [ ...
    "^IXIC"
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
    ];

Name = [ ...
    "Nasdaq Composite"
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
    ];

T = table(Symbol, Name);

end

% -------------------------------------------------------------------------
function T = localNYSE()

Symbol = [ ...
    "^NYA"
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
    ];

Name = [ ...
    "NYSE Composite"
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
    ];

T = table(Symbol, Name);

end

% -------------------------------------------------------------------------
function T = localMilan()

Symbol = [ ...
    "FTSEMIB.MI"
    "ENEL.MI"
    "ENI.MI"
    "ISP.MI"
    "UCG.MI"
    "G.MI"
    "STMMI.MI"
    "TIT.MI"
    "SRG.MI"
    "MB.MI"
    "PST.MI"
    "MONC.MI"
    "LDO.MI"
    "CPR.MI"
    "REC.MI"
    "BZU.MI"
    "BAMI.MI"
    "TEN.MI"
    "HER.MI"
    "A2A.MI"
    "INW.MI"
    "NEXI.MI"
    "RACE.MI"
    "IVG.MI"
    "ERG.MI"
    "ACE.MI"
    "PIRC.MI"
    "DIA.MI"
    "AMP.MI"
    "FBK.MI"
    ];

Name = [ ...
    "FTSE MIB Index"
    "Enel"
    "Eni"
    "Intesa Sanpaolo"
    "UniCredit"
    "Assicurazioni Generali"
    "STMicroelectronics"
    "Telecom Italia"
    "Snam"
    "Mediobanca"
    "Poste Italiane"
    "Moncler"
    "Leonardo"
    "Campari"
    "Recordati"
    "Buzzi"
    "Banco BPM"
    "Tenaris"
    "Hera"
    "A2A"
    "Inwit"
    "Nexi"
    "Ferrari"
    "Iveco Group"
    "ERG"
    "Acea"
    "Pirelli"
    "DiaSorin"
    "Amplifon"
    "FinecoBank"
    ];

T = table(Symbol, Name);

end

% -------------------------------------------------------------------------
function T = localLondon()

Symbol = [ ...
    "^FTSE"
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
    ];

Name = [ ...
    "FTSE 100 Index"
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
    ];

T = table(Symbol, Name);

end

function T = localSP500()

Symbol = [ ...
    "^GSPC"
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
    ];

Name = [ ...
    "S&P 500 Index"
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
    ];

T = table(Symbol, Name);

end