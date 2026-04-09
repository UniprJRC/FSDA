function [out]=getYahoo(ticker, varargin)
%getYahoo downloads financial time series from Yahoo Finance and optionally plots them using a customizable three-panel layout
%
%
%<a href="matlab: docsearchFS('getYahoo')">Link to the help function</a>
%
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
% ticker :      Ticker symbol(s). Character, string or cell array of char.
%               It can be:
%               - a character scalar, for example 'G.MI'
%               - a string scalar or string array, for example "G.MI" or ["G.MI" "ENEL.MI"]
%               - a cell array of character vectors, for example {'G.MI','ENEL.MI'}
%               The function downloads OHLCV data from Yahoo Finance for
%               each supplied ticker and returns the results inside a
%               structure array.
%               Example - 'G.MI'
%               Data Types - char | string | cell
%
% Optional input arguments:
%
%  LastPeriod : Period to analyze. Character or string scalar.
%               Possible values are:
%               '1d' '5d' '1mo' '3mo' '6mo' '1y' '2y' '5y' '10y' 'ytd' 'max'
%               Default is '1y'.
%               Example - 'LastPeriod','6mo'
%               Data Types - char | string
%
%   interval  : Sampling interval requested to Yahoo Finance.
%               Character or string scalar.
%               Possible values are:
%               '1m','2m','5m','15m','30m','60m','90m','1h','1d','5d','1wk','1mo','3mo'
%               Default is '1m'. If the pair LastPeriod/interval is not
%               admissible, the interval is replaced automatically if
%               option autoFixInterval=true.
%               Example - 'interval','1d'
%               Data Types - char | string
%
% widthFactor : Candle width scaling factor. Scalar.
%               It is used only when topPanelMode='candle'.
%               Default is 1.5.
%               Example - 'widthFactor',1.2
%               Data Types - double
%
% layoutHeights : Relative heights of the three panels. Numeric vector.
%               Vector with three positive elements controlling the
%               relative heights of:
%               1) top panel   = price / candles / moving averages;
%               2) middle panel= volume;
%               3) bottom panel= technical indicator.
%               The actual heights are obtained by normalizing the vector
%               so that its elements sum to 1.
%               Default is [1 1 1], which gives the same height to the
%               three panels.
%               For example:
%               - [2 1 1] makes the top panel occupy about half of the
%                 plotting area and the other two panels about one quarter
%                 each.
%               - [3 1 1] gives even more emphasis to the top panel.
%               - [1.5 1 0.7] gives moderate emphasis to the top panel and
%                 reduces the space of the bottom panel.
%               Example - 'layoutHeights',[2 1 1]
%               Data Types - double
%
% removeGaps  : Remove night or weekend gaps in x-axis. Boolean.
%               If true, x-axis is based on the progressive index of the
%               time series.
%               If false, x-axis is based on actual datetime.
%               Default is true.
%               Example - 'removeGaps',false
%               Data Types - logical
%
% breakAtSession : Add separators between sessions. Boolean.
%               Used only when removeGaps=true.
%               Default is true.
%               Example - 'breakAtSession',false
%               Data Types - logical
%
%    nTicks   : Number of ticks shown on x-axis when removeGaps=true.
%               Positive integer scalar.
%               Default is 8.
%               Example - 'nTicks',10
%               Data Types - double
%
% topPanelMode : Type of plot in the first panel. String or char.
%               Possible values are:
%               'candle' = candlestick chart
%               'line'   = close price only
%               'ma'     = close price and moving averages
%               Default is 'candle'.
%               Example - 'topPanelMode','ma'
%               Data Types - char | string
%
% maFastLen   : Fast moving average length. Positive scalar integer.
%               Used only when topPanelMode='ma'.
%               Default is 10.
%               Example - 'maFastLen',5
%               Data Types - double
%
% maMidLen    : Medium moving average length. Positive scalar integer.
%               Used only when topPanelMode='ma'.
%               Default is 30.
%               Example - 'maMidLen',20
%               Data Types - double
%
% maSlowLen   : Slow moving average length. Positive scalar integer.
%               Used only when topPanelMode='ma'.
%               Default is 60. Note that maSlowLen > maMidLen > maFastLen.
%               Example - 'maSlowLen',100
%               Data Types - double
%
% showMACrossovers : Show moving average crossover markers. Boolean.
%               Used only when topPanelMode='ma'.
%               Default is true.
%               Example - 'showMACrossovers',false
%               Data Types - logical
%
% upColor     : Color for up candles (top panel) and bars (mid panel). RGB row vector.
%               Default is [0 0.7 0].
%               Example - 'upColor',[0 0.5 0]
%               Data Types - double
%
% downColor   :  Color for down candles (top panel) and bars (mid panel). RGB row vector.
%               Default is [0.85 0 0].
%               Example - 'downColor',[0.7 0 0]
%               Data Types - double
%
% bottomPanelMode : Type of indicator shown in the third panel. String or char.
%               Possible values are:
%               'rsi'       = Relative Strength Index
%               'stoch'     = Stochastic oscillator
%               'macd'      = MACD
%               'williamsr' = Williams %R
%               'roc'       = Rate of Change
%               Default is 'stoch'.
%               Example - 'bottomPanelMode','macd'
%               Data Types - char | string
%
% topRSI      : Upper RSI threshold. Scalar in in the interval [50 100].
%               Default is 80.
%               Example - 'topRSI',70
%               Data Types - double
%
% lowRSI      : Lower RSI threshold. Scalar in in the interval [0 50].
%               Default is 20.
%               Example - 'lowRSI',30
%               Data Types - double
%
% topStoch    : Upper stochastic threshold. Scalar in in the interval [50 100].
%               Default is 70.
%               Example - 'topStoch',80
%               Data Types - double
%
% lowStoch    : Lower stochastic threshold. Scalar in in the interval [0 50].
%               Default is 30.
%               Example - 'lowStoch',20
%               Data Types - double
%
% topWilliams : Upper Williams %R threshold. Scalar in in the interval [-50 0].
%               Default is -30.
%               Example - 'topWilliams',-20
%               Data Types - double
%
% lowWilliams : Lower Williams %R threshold. Scalar in in the interval [-100 -50].
%               Default is -70.
%               Example - 'lowWilliams',-80
%               Data Types - double
%
% rsiLen      : RSI window length. Scalar. Default is 14.
%               Example - 'rsiLen',10
%               Data Types - double
%
% stochLen    : Stochastic lookback length. Scalar. Default is 14.
%               Example - 'stochLen',10
%               Data Types - double
%
% stochSmooth : Stochastic smoothing length for %%D. Scalar. Default is 3.
%               Example - 'stochSmooth',5
%               Data Types - double
%
% macdFastLen : MACD fast EMA length. Positive scalar integer. Default is 12.
%               Example - 'macdFastLen',8
%               Data Types - double
%
% macdSlowLen : MACD slow EMA length. Positive scalar integer. Default is 26.
%               Example - 'macdSlowLen',17
%               Data Types - double
%
% macdSigLen  : MACD signal EMA length. Positive scalar integer. Default is 9.
%               Example - 'macdSigLen',5
%               Data Types - double
%
% rocLen      : ROC lag length. Positive scalar integer. Default is 12.
%               Example - 'rocLen',6
%               Data Types - double
%
% wrLen       : Williams %%R lookback length. Scalar. Default is 14.
%               Example - 'wrLen',10
%               Data Types - double
%
% showPanelHelp : Display textual explanation in Command Window. Boolean.
%               Default is false.
%               Example - 'showPanelHelp',true
%               Data Types - logical
%
% autoFixInterval : Replace invalid range/interval combinations automatically.
%               Boolean. Default is true.
%               Example - 'autoFixInterval',false
%               Data Types - logical
%
% plots       : Produce plots. Boolean.
%               If true (default) the three-panel figure is created for each
%               ticker. If false only the data are downloaded and stored in
%               the output structure.
%               Example - 'plots',false
%               Data Types - logical
%
% msg         : Display progress messages. Boolean.
%               If true (default), progress messages and warnings are shown.
%               Example - 'msg',false
%               Data Types - logical
%
% Output:
%
% out         : structure array containing the following fields
%
% out.Ticker            = ticker symbol.
% out.LastPeriod        = requested period.
% out.intervalRequested = requested interval after validation.
% out.intervalActual    = actual data granularity returned by Yahoo.
% out.TimeZone          = exchange timezone.
% out.TT                = timetable with variables Open, High, Low, Close, Volume.
% out.Indicators        = structure containing RSI, StochK, StochD, MACD,
%                         MACDSignal, MACDHist, WilliamsR, ROC, maFast,
%                         maMid and maSlow.
% out.Success           = true if download and processing succeeded.
% out.Message           = message describing the result.
% out.class             = 'getYahoo'.
%
% See also: getTickers, getFundamentals, rsindex, candle, movavg
%
% References:
%
%  Damodaran, A. (2012), "Investment Valuation: Tools and Techniques for
%  Determining the Value of Any Asset", 3rd Edition, Wiley.
%
%  Cochrane, J. H. (2023), "Asset Pricing, Revised Edition",
%  Princeton University Press, Princeton.
%
%  Koller, T., Goedhart, M., and Wessels, D. (2020),
%  "Valuation: Measuring and Managing the Value of Companies,
%  7th Edition", Wiley, Hoboken.
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('getYahoo')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Single ticker with default options.
    out=getYahoo('G.MI');
%}

%{
    % Two tickers, no plots.
    out=getYahoo({'G.MI','ENEL.MI'},'plots',false);
%}

%{
    % Daily data for one year with MACD in the bottom panel.
     out=getYahoo("ENEL.MI",'interval','1d','bottomPanelMode','macd');
%}

%{
    % Close price with moving averages in the top panel.
    % Bottom panel hosts the RSI index.
    out=getYahoo("G.MI",'topPanelMode','ma','bottomPanelMode','rsi');
%}

%{
    %% Multiple tickers.
    % ticker passed as a cell array of characters.
    ticker={'G.MI','ENEL.MI','ISP.MI'};
    out = getYahoo(ticker,'plots',false,'msg',false);
%}

%{
    % Auto interval fix.
    out = getYahoo('ENEL.MI','LastPeriod','1y','interval','1m','plots',false);
    
    % Note that with 'autoFixInterval',false the previous examples fails
    % out = getYahoo('ENEL.MI','LastPeriod','1y','interval','1m','autoFixInterval',false,'plots',false);
%}

%{
    % Top panel occupies about half of the available plotting area.
    out=getYahoo('G.MI','layoutHeights',[2 1 1]);
%}

%{
    % Strong emphasis on the price panel with RSI in the bottom panel.
    out=getYahoo('ENEL.MI','topPanelMode','candle', ...
        'bottomPanelMode','rsi','layoutHeights',[3 1 1]);
%}

%{
    % Balanced layout with slightly smaller bottom indicator panel.
    out=getYahoo('G.MI','topPanelMode','ma', ...
        'bottomPanelMode','macd','layoutHeights',[1.5 1 0.7]);
%}

%{
    % Intraday data with real time axis and larger top panel.
    out=getYahoo('G.MI','LastPeriod','5d','interval','15m', ...
        'removeGaps',false,'layoutHeights',[2 1 1]);
%}

%{
    % Monthly data over many years with Williams %%R and custom panel heights.
    out=getYahoo('ENEL.MI','LastPeriod','10y','interval','1mo', ...
        'bottomPanelMode','williamsr','layoutHeights',[2.5 1 1]);
%}

%{
    % Intraday with real time axis.
    out = getYahoo('G.MI','LastPeriod','5d','interval','15m','removeGaps',false);
%}

%{
    % Intraday compressed.
    out = getYahoo('G.MI','LastPeriod','5d','interval','15m','removeGaps',true,'breakAtSession',false);
%}

%{
    % More ticks on the x axis.
    out = getYahoo('G.MI','LastPeriod','3mo','interval','1d','nTicks',12);
%}

%{
    % Line chart in the top panel.
    out = getYahoo('ENEL.MI','topPanelMode','line','bottomPanelMode','rsi');
%}

%{
    %% Line Chart with moving averages in the top panel and crossovers.
    out = getYahoo('ENEL.MI','topPanelMode','ma','maFastLen',5,'maMidLen',20,'maSlowLen',50);
%}

%{
    % Line Chart with moving averages in the top panel but no crossovers.
    out = getYahoo('ENEL.MI','topPanelMode','ma','showMACrossovers',false);
%}

%{
    % Custom colors for the candles and the bars associted with the quantities.
    out = getYahoo('G.MI','upColor',[0 0.4 0.8],'downColor',[0.8 0.4 0]);
%}

%{
    %  Wider candles in the top panel.
    out = getYahoo('G.MI','widthFactor',2);
%}

%{
    % RSI custom in the bottom panel.
    out = getYahoo('G.MI','bottomPanelMode','rsi','rsiLen',10,'topRSI',70,'lowRSI',30);
%}

%{
    % Stochastic custom in the bottom panel.
    out = getYahoo('G.MI','bottomPanelMode','stoch','stochLen',10,'stochSmooth',5);
%}

%{
    % MACD custom in the bottom panel.
    out = getYahoo('ENEL.MI','bottomPanelMode','macd','macdFastLen',8,'macdSlowLen',17);
%}

%{
    % Williams R custom in the bottom panel.
    out = getYahoo('G.MI','bottomPanelMode','williamsr','wrLen',100);
%}

%{
    % ROC custom in the bottom panel.
    out = getYahoo('G.MI','bottomPanelMode','roc','rocLen',6);
%}

%{
    % Long horizon of 10 years interval 1mo.
    out = getYahoo('ENEL.MI','LastPeriod','10y','interval','1mo');
%}

%{
    % Options LastPeriod, with ytd.
    out = getYahoo('ENEL.MI','LastPeriod','ytd','interval','1d');
%}

%{
    %% Max history.
    out = getYahoo('G.MI','LastPeriod','max','interval','1wk','plots',false);
%}

%{
    %% With help.
    out = getYahoo('G.MI','showPanelHelp',true);
%}

%{
    % Inspect output structure array.
    tickers=["G.MI","ENEL.MI"];
    out = getYahoo(tickers,'plots',false);
    disp(tickers(1))
    disp(out(1))
    disp(tickers(2))
    disp(out(2))
%}

%{
    % Compare indicators in the bottom panel.
    getYahoo('G.MI','bottomPanelMode','rsi');
    getYahoo('G.MI','bottomPanelMode','macd');
%}

%{
    % Intraday 1d with 1m.
    out = getYahoo('ENEL.MI','LastPeriod','1d','interval','1m');
%}

%{
    % Full options stress test.
    out = getYahoo('G.MI','LastPeriod','3mo','interval','1d','topPanelMode','ma',...
        'maFastLen',7,'maMidLen',21,'maSlowLen',50,'bottomPanelMode','rsi');
%}

%{
    % Use output.
    ticker='G.MI';
    out = getYahoo(ticker,'plots',false);
    TT = out.TT;
    plot(TT.t,TT.Close);
    title(ticker)
%}


%{
    % Combined use of getTickers, getYahoo and getFundamentals.
    T = getTickers('market','London','Source','dynamic','nStocks',15);
    disp(T)

    % Retrieve price data
    out = getYahoo(T.ticker(2:end));

    % Retrieve valuation metrics
    F = getFundamentals(T.ticker(2:end),'Fields','valuation');
    disp(F)
%}

%% Beginning of code

if nargin<1
    error('FSDA:getYahoo:MissingInput','Ticker symbol is missing.');
end

if coder.target('MATLAB')
    options=struct( ...
        'LastPeriod','1y', ...
        'interval','1m', ...
        'widthFactor',1.5, ...
        'removeGaps',true, ...
        'breakAtSession',true, ...
        'nTicks',8, ...
        'topPanelMode','candle', ...
        'maFastLen',10, ...
        'maMidLen',30, ...
        'maSlowLen',60, ...
        'showMACrossovers',true, ...
        'upColor',[0 0.7 0], ...
        'downColor',[0.85 0 0], ...
        'bottomPanelMode','stoch', ...
        'topRSI',80, ...
        'lowRSI',20, ...
        'topStoch',70, ...
        'lowStoch',30, ...
        'topWilliams',-30, ...
        'lowWilliams',-70, ...
        'rsiLen',14, ...
        'stochLen',14, ...
        'stochSmooth',3, ...
        'macdFastLen',12, ...
        'macdSlowLen',26, ...
        'macdSigLen',9, ...
        'rocLen',12, ...
        'wrLen',14, ...
        'showPanelHelp',false, ...
        'autoFixInterval',true, ...
        'plots',true, ...
        'layoutHeights',[1 1 1], ...
        'msg',true);

    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:getYahoo:WrongInputOpt', ...
                'Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        aux.chkoptions(options,UserOptions);
    end
end

if nargin>1
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

LastPeriod       = char(string(options.LastPeriod));
interval         = char(string(options.interval));
widthFactor      = options.widthFactor;
removeGaps       = options.removeGaps;
breakAtSession   = options.breakAtSession;
nTicks           = options.nTicks;
topPanelMode     = char(string(options.topPanelMode));
maFastLen        = options.maFastLen;
maMidLen         = options.maMidLen;
maSlowLen        = options.maSlowLen;
showMACrossovers = options.showMACrossovers;
upColor          = options.upColor;
downColor        = options.downColor;
bottomPanelMode  = char(string(options.bottomPanelMode));
topRSI           = options.topRSI;
lowRSI           = options.lowRSI;
topStoch         = options.topStoch;
lowStoch         = options.lowStoch;
topWilliams      = options.topWilliams;
lowWilliams      = options.lowWilliams;
rsiLen           = options.rsiLen;
stochLen         = options.stochLen;
stochSmooth      = options.stochSmooth;
macdFastLen      = options.macdFastLen;
macdSlowLen      = options.macdSlowLen;
macdSigLen       = options.macdSigLen;
rocLen           = options.rocLen;
wrLen            = options.wrLen;
showPanelHelp    = options.showPanelHelp;
autoFixInterval  = options.autoFixInterval;
doplot           = options.plots;
msg              = options.msg;
layoutHeights    = options.layoutHeights;

if isempty(interval)
    switch LastPeriod
        case '1d'
            interval = '1m';
        case '5d'
            interval = '5m';
        case {'1mo','3mo'}
            interval = '30m';
        otherwise
            interval = '1d';
    end
end

% Normalize ticker input to string array
if ischar(ticker)
    tickerList = string({ticker});
elseif isstring(ticker)
    tickerList = ticker(:);
elseif iscell(ticker) && all(cellfun(@ischar,ticker))
    tickerList = string(ticker);
else
    error('FSDA:getYahoo:WngInp', ...
        'ticker must be a char, a cell array of char, or a string array.');
end

% Validate range/interval pair
[intervalValidated, isValidPair] = validateYahooRangeInterval(LastPeriod, interval, autoFixInterval);

if ~isValidPair && msg
    fprintf('Requested pair range=%s, interval=%s is not supported.\n', LastPeriod, interval);
    fprintf('Using interval=%s instead.\n', intervalValidated);
end

hf = ["User-Agent", "Mozilla/5.0"; ...
    "Accept", "application/json,text/plain,*/*"];
opts = weboptions('ContentType','json','Timeout',30,'HeaderFields',hf);

neededFields = {'open','high','low','close','volume'};

nout = numel(tickerList);
out = repmat(struct( ...
    'Ticker',"", ...
    'LastPeriod',"", ...
    'intervalRequested',"", ...
    'intervalActual',"", ...
    'TimeZone',"", ...
    'TT',timetable, ...
    'Indicators',struct, ...
    'Success',false, ...
    'Message',"", ...
    'class','getYahoo'), nout,1);

for ii=1:nout
    tickerNow = tickerList(ii);

    out(ii).Ticker            = tickerNow;
    out(ii).LastPeriod        = string(LastPeriod);
    out(ii).intervalRequested = string(intervalValidated);
    out(ii).class             = 'getYahoo';

    if msg
        fprintf('\nProcessing ticker %d of %d: %s\n', ii, nout, tickerNow);
    end

    url = "https://query1.finance.yahoo.com/v8/finance/chart/" + tickerNow;

    try
        raw = webread(url, 'interval', intervalValidated, 'range', LastPeriod, opts);
    catch ME
        out(ii).Message = string(ME.message);
        warning('FSDA:getYahoo:RequestFailed', ...
            'Yahoo request failed for ticker=%s, range=%s, interval=%s.\nOriginal error: %s', ...
            tickerNow, LastPeriod, intervalValidated, ME.message);
        continue
    end

    if ~isfield(raw,'chart') || ~isfield(raw.chart,'result')
        out(ii).Message = "Yahoo response does not contain chart.result.";
        warning('FSDA:getYahoo:MalformedReply', ...
            'Yahoo returned a malformed reply for ticker=%s.', tickerNow);
        continue
    end

    r = raw.chart.result;
    if isempty(r)
        out(ii).Message = "Yahoo returned an empty result.";
        warning('FSDA:getYahoo:EmptyResult', ...
            'Yahoo returned an empty result for ticker=%s.', tickerNow);
        continue
    elseif iscell(r)
        r = r{1};
    else
        r = r(1);
    end

    if ~isfield(r,'meta') || ~isfield(r.meta,'dataGranularity') || ~isfield(r.meta,'exchangeTimezoneName')
        out(ii).Message = "Yahoo returned incomplete metadata.";
        warning('FSDA:getYahoo:MissingMeta', ...
            'Yahoo returned incomplete metadata for ticker=%s.', tickerNow);
        continue
    end

    granularity = string(r.meta.dataGranularity);
    intervalThis = granularity;
    out(ii).intervalActual = intervalThis;
    out(ii).TimeZone = string(r.meta.exchangeTimezoneName);

    if msg && string(intervalValidated) ~= granularity
        fprintf('Requested time interval is: %s for ticker %s\n', intervalValidated, tickerNow);
        fprintf('Downloaded time interval is: %s for ticker %s\n', intervalThis, tickerNow);
    end

    tz = r.meta.exchangeTimezoneName;

    if ~isfield(r,'timestamp') || isempty(r.timestamp)
        out(ii).Message = "Yahoo returned empty timestamps.";
        warning('FSDA:getYahoo:EmptyTimestamp', ...
            'Yahoo returned empty timestamps for ticker=%s.', tickerNow);
        continue
    end

    t = datetime(r.timestamp, 'ConvertFrom','posixtime', 'TimeZone', tz);

    if ~isfield(r,'indicators') || ~isfield(r.indicators,'quote')
        out(ii).Message = "Yahoo returned no quote field.";
        warning('FSDA:getYahoo:MissingQuote', ...
            'Yahoo returned no quote field for ticker=%s.', tickerNow);
        continue
    end

    q = r.indicators.quote;
    if isempty(q)
        out(ii).Message = "Yahoo returned empty quote data.";
        warning('FSDA:getYahoo:EmptyQuote', ...
            'Yahoo returned empty quote data for ticker=%s.', tickerNow);
        continue
    elseif iscell(q)
        q = q{1};
    else
        q = q(1);
    end

    if ~all(isfield(q,neededFields))
        out(ii).Message = "Yahoo returned incomplete OHLCV fields.";
        warning('FSDA:getYahoo:MissingOHLCV', ...
            'Yahoo returned incomplete OHLCV fields for ticker=%s.', tickerNow);
        continue
    end

    n = min([numel(t), numel(q.open), numel(q.high), numel(q.low), ...
        numel(q.close), numel(q.volume)]);

    if n==0
        out(ii).Message = "Yahoo returned empty OHLCV series.";
        warning('FSDA:getYahoo:EmptySeries', ...
            'Yahoo returned empty OHLCV series for ticker=%s.', tickerNow);
        continue
    end

    t  = t(1:n);        t  = t(:);
    op = q.open(1:n);   op = op(:);
    hi = q.high(1:n);   hi = hi(:);
    lo = q.low(1:n);    lo = lo(:);
    cl = q.close(1:n);  cl = cl(:);
    vo = q.volume(1:n); vo = vo(:);

    TT = timetable(t, op, hi, lo, cl, vo, ...
        'VariableNames', {'Open','High','Low','Close','Volume'});

    TT = rmmissing(TT);

    if isempty(TT)
        out(ii).Message = "No valid observations available after removing missing values.";
        warning('FSDA:getYahoo:EmptyData', ...
            'No valid observations available for ticker=%s after removing missing values.', tickerNow);
        continue
    end

    % Indicators
    if height(TT) > 5
        RSI = rsindex(TT,'WindowSize',rsiLen);
        rsiVals = RSI.RelativeStrengthIndex;
    else
        rsiVals = 50 * ones(height(TT),1);
    end

    LL = movmin(TT.Low,[stochLen-1 0],'omitnan');
    HH = movmax(TT.High,[stochLen-1 0],'omitnan');
    denStoch = HH - LL;
    stochK = nan(height(TT),1);
    ok = denStoch ~= 0;
    stochK(ok) = 100 * (TT.Close(ok) - LL(ok)) ./ denStoch(ok);
    stochD = movmean(stochK, stochSmooth, 'omitnan');

    emaFast = movavg(TT.Close, 'exponential', macdFastLen);
    emaSlow = movavg(TT.Close, 'exponential', macdSlowLen);
    macdLine   = emaFast - emaSlow;
    macdSignal = movavg(macdLine, 'exponential', macdSigLen);
    macdHist   = macdLine - macdSignal;

    wrHH = movmax(TT.High,[wrLen-1 0],'omitnan');
    wrLL = movmin(TT.Low,[wrLen-1 0],'omitnan');
    denWR = wrHH - wrLL;
    williamsR = nan(height(TT),1);
    okWR = denWR ~= 0;
    williamsR(okWR) = -100 * (wrHH(okWR) - TT.Close(okWR)) ./ denWR(okWR);

    rocVals = nan(height(TT),1);
    if height(TT) > rocLen
        closeLag = [nan(rocLen,1); TT.Close(1:end-rocLen)];
        okROC = closeLag ~= 0 & ~isnan(closeLag);
        rocVals(okROC) = 100 * (TT.Close(okROC) - closeLag(okROC)) ./ closeLag(okROC);
    end

    maFast = movmean(TT.Close, maFastLen, 'omitnan');
    maMid  = movmean(TT.Close, maMidLen,  'omitnan');
    maSlow = movmean(TT.Close, maSlowLen, 'omitnan');

    out(ii).TT = TT;
    out(ii).Indicators = struct( ...
        'RSI',rsiVals, ...
        'StochK',stochK, ...
        'StochD',stochD, ...
        'MACD',macdLine, ...
        'MACDSignal',macdSignal, ...
        'MACDHist',macdHist, ...
        'WilliamsR',williamsR, ...
        'ROC',rocVals, ...
        'maFast',maFast, ...
        'maMid',maMid, ...
        'maSlow',maSlow);
    out(ii).Success = true;
    out(ii).Message = "OK";

    if doplot
        localPlotYahoo(TT, tickerNow, LastPeriod, intervalThis, ...
            widthFactor, removeGaps, breakAtSession, nTicks, ...
            topPanelMode, maFastLen, maMidLen, maSlowLen, showMACrossovers, ...
            upColor, downColor, ...
            bottomPanelMode, topRSI, lowRSI, topStoch, lowStoch, ...
            topWilliams, lowWilliams, ...
            rsiVals, stochK, stochD, macdLine, macdSignal, macdHist, ...
            williamsR, rocVals, maFast, maMid, maSlow,layoutHeights)
    end
end

if showPanelHelp
    showPanelExplanationCommand(topPanelMode, bottomPanelMode, ...
        maFastLen, maMidLen, maSlowLen, ...
        rsiLen, stochLen, stochSmooth, macdFastLen, macdSlowLen, macdSigLen, rocLen, wrLen, ...
        topRSI, lowRSI, topStoch, lowStoch, topWilliams, lowWilliams);
end

end
function localPlotYahoo(TT, tickerNow, LastPeriod, intervalThis, ...
    widthFactor, removeGaps, breakAtSession, nTicks, ...
    topPanelMode, maFastLen, maMidLen, maSlowLen, showMACrossovers, ...
    upColor, downColor, ...
    bottomPanelMode, topRSI, lowRSI, topStoch, lowStoch, ...
    topWilliams, lowWilliams, ...
    rsiVals, stochK, stochD, macdLine, macdSignal, macdHist, ...
    williamsR, rocVals, maFast, maMid, maSlow, layoutHeights)

layoutHeights = layoutHeights(:)';
if numel(layoutHeights) ~= 3 || any(~isfinite(layoutHeights)) || any(layoutHeights <= 0)
    error('FSDA:getYahoo:WngInp', ...
        'layoutHeights must be a numeric vector with 3 positive elements.');
end

layoutHeights = layoutHeights / sum(layoutHeights);

% X coordinates
if removeGaps
    x = (1:height(TT))';
else
    x = TT.t;
end

% Break indicator lines when removeGaps = false
if ~removeGaps
    dt = diff(TT.t);

    if strcmp(intervalThis,'1d')
        gapIdx = [false; hours(dt) > 40];
    else
        gapIdx = [false; hours(dt) > 10];
    end

    stochK(gapIdx)     = NaN;
    stochD(gapIdx)     = NaN;
    rsiVals(gapIdx)    = NaN;
    macdLine(gapIdx)   = NaN;
    macdSignal(gapIdx) = NaN;
    macdHist(gapIdx)   = NaN;
    williamsR(gapIdx)  = NaN;
    rocVals(gapIdx)    = NaN;
end

% Extra info for datatips
dirTxt = strings(height(TT),1);
dirTxt(TT.Close >= TT.Open) = "Up candle";
dirTxt(TT.Close <  TT.Open) = "Down candle";

retPct = nan(height(TT),1);
prevClose = TT.Close(1:end-1);
dClose = diff(TT.Close);
okRet = prevClose ~= 0 & ~isnan(prevClose);
idxRet = find(okRet) + 1;
retPct(idxRet) = 100 * dClose(okRet) ./ prevClose(okRet);

% Figure and manual layout
figure('Color','w');

leftMargin   = 0.07;
rightMargin  = 0.03;
bottomMargin = 0.08;
topMargin    = 0.08;
vGap         = 0.045;

availW = 1 - leftMargin - rightMargin;
availH = 1 - topMargin - bottomMargin - 2*vGap;

h1 = availH * layoutHeights(1);
h2 = availH * layoutHeights(2);
h3 = availH * layoutHeights(3);

y3 = bottomMargin;
y2 = y3 + h3 + vGap;
y1 = y2 + h2 + vGap;

ax1 = axes('Position',[leftMargin y1 availW h1]);
hold(ax1,'on');
grid(ax1,'on');
ylabel(ax1,'Price');

%% Top panel
switch lower(topPanelMode)
   case 'bollinger'
    % Parametri Bollinger
    windowLen = 20;
    nSigma = 2;

    % Assicurati di usare le stesse coordinate x viste dal resto della funzione
    if removeGaps
        xPlot = (1:height(TT))';
    else
        xPlot = TT.t;
    end

    priceAll = TT.Close;      % mantiene NaN dove presenti
    validIdx = ~isnan(priceAll);

    if nnz(validIdx) < 2
        warning('FSDA:getYahoo:BollingerNoData','Troppi pochi dati validi per Bollinger per %s.', tickerNow);
        plot(ax1, xPlot, priceAll, '-', 'Color', [0.2 0.2 0.2]);
    else
    end

    % Calcola SMA e std in modo allineato alla serie completa:
    % movmean/movstd con finestra causal (windowLen) e 'omitnan'
    if numel(priceAll) < windowLen
        smaAll = movmean(priceAll, [numel(priceAll)-1 0], 'omitnan');
        stdAll = movstd(priceAll, [numel(priceAll)-1 0], 'omitnan');
    else
        smaAll = movmean(priceAll, [windowLen-1 0], 'omitnan');
        stdAll = movstd(priceAll, [windowLen-1 0], 'omitnan');
    end
    upperAll = smaAll + nSigma .* stdAll;
    lowerAll = smaAll - nSigma .* stdAll;

    % Traccia prezzo con colori up/down: determina up/down su dati validi
    dPrice = [NaN; diff(priceAll)];
    isUp = dPrice >= 0;

    % Per plotting a segmenti continui, disegniamo dove price non NaN
    plot(ax1, xPlot(validIdx & isUp(validIdx)), priceAll(validIdx & isUp(validIdx)), '-', ...
         'Color', upColor, 'LineWidth', 1.2);
    plot(ax1, xPlot(validIdx & ~isUp(validIdx)), priceAll(validIdx & ~isUp(validIdx)), '-', ...
         'Color', downColor, 'LineWidth', 1.2);

    % SMA e bande (mostriamo anche punti NaN come gap)
    plot(ax1, xPlot, smaAll, '-','Color',[0 0.4470 0.7410],'LineWidth',1);
    hUp = plot(ax1, xPlot, upperAll, '--','Color',[0.3 0.3 0.3],'LineWidth',0.9);
    hLo = plot(ax1, xPlot, lowerAll, '--','Color',[0.3 0.3 0.3],'LineWidth',0.9);

    % Riempimento tra bande (funziona sia con datetime sia con indici,
    % rimuovendo punti NaN per patch)
    okPatch = ~isnan(upperAll) & ~isnan(lowerAll) & ~isnan(xPlot);
    if nnz(okPatch) >= 2
        xPatch = [xPlot(okPatch); flipud(xPlot(okPatch))];
        yPatch = [upperAll(okPatch); flipud(lowerAll(okPatch))];
        patch('XData', xPatch, 'YData', yPatch, 'FaceColor',[0.6 0.6 0.6], ...
              'FaceAlpha',0.08, 'EdgeColor','none', 'Parent', ax1);
    end

    % Evidenzia ultimo punto valido
    lastIdx = find(validIdx, 1, 'last');
    if ~isempty(lastIdx)
        plot(ax1, xPlot(lastIdx), priceAll(lastIdx), 'o', ...
             'MarkerFaceColor',[0 0.4470 0.7410], 'MarkerEdgeColor','k');
    end

    % Limiti: usa tutti i dati validi (prezzo e bande)
    yrCandidates = [priceAll(validIdx); upperAll(okPatch); lowerAll(okPatch)];
    if ~isempty(yrCandidates)
        yr = [min(yrCandidates) max(yrCandidates)];
        pad = 0.06 * (yr(2)-yr(1));
        if pad == 0
            pad = max(abs(yr(2)),1)*0.01;
        end
        ylim(ax1, [yr(1)-pad, yr(2)+pad]);
    end

    legend(ax1, {'Price up','Price down','SMA','Upper','Lower'}, 'Location','best');

    case 'candle'
        if removeGaps
            Data = [TT.Open TT.High TT.Low TT.Close];
            h = candle(Data,'b');
        else
            h = candle(TT);
            Data = TT{:,1:4};
        end

        patches = findobj(ax1,'Type','Patch');
        lines   = findobj(ax1,'Type','Line');

        opens  = Data(:,1);
        closes = Data(:,4);
        isUp   = closes >= opens;

        patches = flip(patches);

        for k = 1:min(numel(patches), numel(isUp))
            if isUp(k)
                set(patches(k),'FaceColor',upColor,'EdgeColor',upColor);
            else
                set(patches(k),'FaceColor',downColor,'EdgeColor',downColor);
            end
        end

        set(lines,'Color',[0 0 0]);

        patchesH = findobj(h,'Type','Patch');
        for k = 1:numel(patchesH)
            xd = get(patchesH(k),'XData');
            xc = mean(xd);
            xnew = xc + widthFactor*(xd - xc);
            set(patchesH(k),'XData',xnew);
        end

    case 'line'
        plot(ax1,x,TT.Close,'k','LineWidth',1.2);

    case 'ma'
        plot(ax1,x,TT.Close,'Color',[0.2 0.2 0.2],'LineWidth',1.0);
        plot(ax1,x,maFast,'b-','LineWidth',1.5);
        plot(ax1,x,maMid,'r--','LineWidth',1.3);
        plot(ax1,x,maSlow,'m:','LineWidth',1.2);

        legend(ax1, ...
            {'Close',sprintf('MA%d',maFastLen),sprintf('MA%d',maMidLen),sprintf('MA%d',maSlowLen)}, ...
            'Location','best');

        if showMACrossovers
            dFM = maFast - maMid;
            bullFM = find(dFM(1:end-1) <= 0 & dFM(2:end) > 0) + 1;
            bearFM = find(dFM(1:end-1) >= 0 & dFM(2:end) < 0) + 1;

            plot(ax1,x(bullFM),maFast(bullFM),'^', ...
                'MarkerSize',7, ...
                'MarkerFaceColor',[0 0.6 0], ...
                'MarkerEdgeColor','k', ...
                'LineStyle','none', ...
                'HandleVisibility','off');

            plot(ax1,x(bearFM),maFast(bearFM),'v', ...
                'MarkerSize',7, ...
                'MarkerFaceColor',[0.85 0 0], ...
                'MarkerEdgeColor','k', ...
                'LineStyle','none', ...
                'HandleVisibility','off');
        end

    otherwise
        error('FSDA:getYahoo:WngInp','Unknown topPanelMode.');
end

ylim(ax1,'tight');

%% Middle panel: volume
ax2 = axes('Position',[leftMargin y2 availW h2]);
hold(ax2,'on');
grid(ax2,'on');
ylabel(ax2,'Volume');

up = TT.Close >= TT.Open;
volColors = zeros(height(TT),3);
volColors(up,:)  = repmat(upColor,   sum(up), 1);
volColors(~up,:) = repmat(downColor, sum(~up), 1);

b = bar(ax2, x, TT.Volume, 'FaceColor','flat','EdgeColor','none');
b.CData = volColors;

%% Bottom panel
ax3 = axes('Position',[leftMargin y3 availW h3]);
hold(ax3,'on');
grid(ax3,'on');
xlabel(ax3,'Time');

switch lower(bottomPanelMode)
    case 'rsi'
        xPatch = [x(1) x(end) x(end) x(1)];
        patch(ax3,xPatch,[topRSI topRSI 100 100],[1 0.9 0.9], ...
            'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');
        patch(ax3,xPatch,[0 0 lowRSI lowRSI],[0.9 0.95 1], ...
            'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');

        hInd = plot(ax3,x,rsiVals,'k','LineWidth',1.5);

        yline(ax3,topRSI,'--r',num2str(topRSI),'LineWidth',1,'HandleVisibility','off');
        yline(ax3,lowRSI,'--b',num2str(lowRSI),'LineWidth',1,'HandleVisibility','off');
        yline(ax3,50,':','50','Color',[0.5 0.5 0.5],'HandleVisibility','off');

        ylabel(ax3,'RSI');
        bottomTitle = "RSI";
        bottomShort = "RSI";
        bottomFixedLims = [0 100];

    case 'stoch'
        xPatch = [x(1) x(end) x(end) x(1)];
        patch(ax3,xPatch,[topStoch topStoch 100 100],[1 0.9 0.9], ...
            'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');
        patch(ax3,xPatch,[0 0 lowStoch lowStoch],[0.9 0.95 1], ...
            'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');

        hK = plot(ax3,x,stochK,'k','LineWidth',1.4);
        hD = plot(ax3,x,stochD,'--','Color',[0.4 0.4 0.4],'LineWidth',1.2);

        yline(ax3,topStoch,'--r',num2str(topStoch),'LineWidth',1,'HandleVisibility','off');
        yline(ax3,lowStoch,'--b',num2str(lowStoch),'LineWidth',1,'HandleVisibility','off');
        yline(ax3,50,':','50','Color',[0.5 0.5 0.5],'HandleVisibility','off');

        ylabel(ax3,'Stochastic (%K, %D)');
        legend(ax3,[hK hD],{'%K','%D'},'Location','best');

        hInd = hK;
        bottomTitle = "Stochastic";
        bottomShort = "%K";
        bottomFixedLims = [0 100];

    case 'macd'
        hBar  = bar(ax3,x,macdHist,'FaceColor',[0.7 0.7 0.7], ...
            'EdgeColor','none','FaceAlpha',0.7);
        hMacd = plot(ax3,x,macdLine,'k','LineWidth',1.4);
        hSig  = plot(ax3,x,macdSignal,'--','Color',[0.5 0 0.8],'LineWidth',1.2);

        yline(ax3,0,':','0','Color',[0.5 0.5 0.5],'HandleVisibility','off');

        ylabel(ax3,'MACD');
        legend(ax3,[hBar hMacd hSig],{'Hist','MACD','Signal'},'Location','best');

        hInd = hMacd;
        bottomTitle = "MACD";
        bottomShort = "MACD";
        bottomFixedLims = [];

    case 'williamsr'
        xPatch = [x(1) x(end) x(end) x(1)];
        patch(ax3,xPatch,[topWilliams topWilliams 0 0],[1 0.9 0.9], ...
            'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');
        patch(ax3,xPatch,[-100 -100 lowWilliams lowWilliams],[0.9 0.95 1], ...
            'EdgeColor','none','FaceAlpha',0.35,'HandleVisibility','off');

        hInd = plot(ax3,x,williamsR,'k','LineWidth',1.5);

        yline(ax3,topWilliams,'--r',num2str(topWilliams),'LineWidth',1,'HandleVisibility','off');
        yline(ax3,lowWilliams,'--b',num2str(lowWilliams),'LineWidth',1,'HandleVisibility','off');
        yline(ax3,-50,':','-50','Color',[0.5 0.5 0.5],'HandleVisibility','off');

        ylabel(ax3,'Williams %R');
        bottomTitle = "Williams %R";
        bottomShort = "W%R";
        bottomFixedLims = [-100 0];

    case 'roc'
        hInd = plot(ax3,x,rocVals,'k','LineWidth',1.5);

        yline(ax3,0,':','0','Color',[0.5 0.5 0.5],'HandleVisibility','off');

        ylabel(ax3,'ROC (%)');
        bottomTitle = "Rate of Change (ROC)";
        bottomShort = "ROC";
        bottomFixedLims = [];

    otherwise
        error('FSDA:getYahoo:WngInp','Unknown bottomPanelMode.');
end

if ~isempty(bottomFixedLims)
    ylim(ax3,bottomFixedLims);
end

sgtitle(sprintf('Price, Volume, and %s for %s (%s, %s)', ...
    bottomTitle, tickerNow, LastPeriod, intervalThis));

%% Datatips
dcm1 = datacursormode(gcf);
set(dcm1,'Enable','on', ...
    'UpdateFcn', @(~,evt) myPriceTip(evt, TT, removeGaps, dirTxt, retPct, ...
    topPanelMode, maFastLen, maMidLen, maSlowLen, maFast, maMid, maSlow, ...
    bottomPanelMode, rsiVals, stochK, stochD, macdLine, macdSignal, macdHist, ...
    williamsR, rocVals));

b.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Date', string(TT.t));
volStr = arrayfun(@fmtVolume, TT.Volume, 'UniformOutput', false);
b.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Volume', volStr);
b.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Close', TT.Close);
addBottomIndicatorRows(b, bottomPanelMode, rsiVals, stochK, stochD, ...
    macdLine, macdSignal, macdHist, williamsR, rocVals);
b.DataTipTemplate.DataTipRows(1).Label = 'X';

hInd.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Date', string(TT.t));
hInd.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Close', TT.Close);
addBottomIndicatorRows(hInd, bottomPanelMode, rsiVals, stochK, stochD, ...
    macdLine, macdSignal, macdHist, williamsR, rocVals);
hInd.DataTipTemplate.DataTipRows(1).Label = 'X';

%% Link axes
linkaxes([ax1 ax2 ax3],'x');

%% X ticks and separators
if removeGaps
    tickPos = round(linspace(1, height(TT), nTicks));
    tickPos = unique(tickPos);

    if ismember(string(intervalThis), ["1m","2m","5m","15m","30m","60m","90m","1h"])
        tickLbl = string(TT.t(tickPos), 'dd-MMM HH:mm');
    else
        tickLbl = string(TT.t(tickPos), 'dd-MMM-yyyy');
    end

    if breakAtSession
        ttLocal = TT.t;
        intervalStr = string(intervalThis);
        idxSep = [];

        nYears = numel(unique(year(ttLocal)));

        if nYears > 1
            yy = year(ttLocal);
            yearChange = [false; diff(yy) ~= 0];
            idxSep = find(yearChange) + 0.5;

        elseif ismember(intervalStr, ["1m","2m","5m","15m","30m","60m","90m","1h"])
            sessionChange = [false; days(diff(ttLocal)) >= 0.5];
            idxSep = find(sessionChange) + 0.5;

        elseif strcmp(intervalStr,'1d')
            ym = year(ttLocal)*100 + month(ttLocal);
            monthChange = [false; diff(ym) ~= 0];
            idxSep = find(monthChange) + 0.5;
        end

        for k = 1:numel(idxSep)
            xline(ax1, idxSep(k), ':', 'Color',[0.75 0.75 0.75], ...
                'LineWidth',2,'HandleVisibility','off');
            xline(ax2, idxSep(k), ':', 'Color',[0.75 0.75 0.75], ...
                'LineWidth',2,'HandleVisibility','off');
            xline(ax3, idxSep(k), ':', 'Color',[0.75 0.75 0.75], ...
                'LineWidth',2,'HandleVisibility','off');
        end
    end

    set([ax1 ax2 ax3],'XTick',tickPos);
    set(ax1,'XTickLabel',[]);
    set(ax2,'XTickLabel',[]);
    set(ax3,'XTickLabel',tickLbl);
else
    ax1.XTickLabel = [];
    ax2.XTickLabel = [];
end

xlim(ax1,'tight');

%% Last labels
switch lower(bottomPanelMode)
    case 'rsi'
        txtVal = rsiVals(end);
    case 'stoch'
        txtVal = stochK(end);
    case 'macd'
        txtVal = macdLine(end);
    case 'williamsr'
        txtVal = williamsR(end);
    case 'roc'
        txtVal = rocVals(end);
end

if ~isnan(txtVal)
    text(ax3, x(end), txtVal, sprintf('  %s\n %.1f', bottomShort, txtVal), ...
        'VerticalAlignment','middle','FontWeight','bold');
end

text(ax1, x(end), TT.Close(end), sprintf('  Last\n %.3f', TT.Close(end)), ...
    'VerticalAlignment','middle','FontWeight','bold');

%% Dynamic rescaling
zh = zoom(gcf);
set(zh,'ActionPostCallback', @(~,evd) rescaleVisibleY(evd, ax1, ax2, ax3, x, TT, ...
    topPanelMode, maFastLen, maMidLen, maSlowLen, ...
    bottomPanelMode, bottomFixedLims, macdLine, macdSignal, macdHist, rocVals));

ph = pan(gcf);
set(ph,'ActionPostCallback', @(~,evd) rescaleVisibleY(evd, ax1, ax2, ax3, x, TT, ...
    topPanelMode, maFastLen, maMidLen, maSlowLen, ...
    bottomPanelMode, bottomFixedLims, macdLine, macdSignal, macdHist, rocVals));

rescaleVisibleY([], ax1, ax2, ax3, x, TT, ...
    topPanelMode, maFastLen, maMidLen, maSlowLen, ...
    bottomPanelMode, bottomFixedLims, macdLine, macdSignal, macdHist, rocVals);

end

function txt = myPriceTip(evt, TT, removeGaps, dirTxt, retPct, ...
    topPanelMode, maFastLen, maMidLen, maSlowLen, maFast, maMid, maSlow, ...
    bottomPanelMode, rsiVals, stochK, stochD, macdLine, macdSignal, macdHist, williamsR, rocVals)

pos = evt.Position;

if removeGaps
    idx = round(pos(1));
    idx = max(1, min(height(TT), idx));
else
    [~, idx] = min(abs(TT.t - pos(1)));
end

switch lower(topPanelMode)
    case 'ma'
        topExtraLines = { ...
            [sprintf('MA%d: ', maFastLen) num2str(maFast(idx), '%.3f')]; ...
            [sprintf('MA%d: ', maMidLen)  num2str(maMid(idx),  '%.3f')]; ...
            [sprintf('MA%d: ', maSlowLen) num2str(maSlow(idx), '%.3f')]};
    otherwise
        topExtraLines = {};
end

switch lower(bottomPanelMode)
    case 'rsi'
        indLines = {['RSI: ' num2str(rsiVals(idx), '%.2f')]};
    case 'stoch'
        indLines = {['Stoch %K: ' num2str(stochK(idx), '%.2f')]; ...
            ['Stoch %D: ' num2str(stochD(idx), '%.2f')]};
    case 'macd'
        indLines = {['MACD: ' num2str(macdLine(idx), '%.4f')]; ...
            ['Signal: ' num2str(macdSignal(idx), '%.4f')]; ...
            ['Hist: ' num2str(macdHist(idx), '%.4f')]};
    case 'williamsr'
        indLines = {['Williams %R: ' num2str(williamsR(idx), '%.2f')]};
    case 'roc'
        indLines = {['ROC: ' num2str(rocVals(idx), '%.2f') '%']};
    otherwise
        indLines = {};
end

txt = [ ...
    {['Date: ' char(datetime(TT.t(idx), 'Format','dd-MMM-yyyy HH:mm'))]
    ['Open: '   num2str(TT.Open(idx),  '%.3f')]
    ['High: '   num2str(TT.High(idx),  '%.3f')]
    ['Low: '    num2str(TT.Low(idx),   '%.3f')]
    ['Close: '  num2str(TT.Close(idx), '%.3f')]}
    topExtraLines
    {['Volume: ' fmtVolume(TT.Volume(idx))]}
    indLines
    {['Type: '   char(dirTxt(idx))]
    ['Return: ' num2str(retPct(idx), '%.2f') '%']}];
end

function rescaleVisibleY(~, ax1, ax2, ax3, x, TT, ...
    topPanelMode, maFastLen, maMidLen, maSlowLen, ...
    bottomPanelMode, bottomFixedLims, macdLine, macdSignal, macdHist, rocVals)

xlimNow = xlim(ax1);
idx = x >= xlimNow(1) & x <= xlimNow(2);

if ~any(idx)
    return
end

switch lower(topPanelMode)
    case 'candle'
        yTopMin = min(TT.Low(idx), [], 'omitnan');
        yTopMax = max(TT.High(idx), [], 'omitnan');
    case 'line'
        yTopMin = min(TT.Close(idx), [], 'omitnan');
        yTopMax = max(TT.Close(idx), [], 'omitnan');
    case 'ma'
        maFast = movmean(TT.Close, maFastLen, 'omitnan');
        maMid  = movmean(TT.Close, maMidLen,  'omitnan');
        maSlow = movmean(TT.Close, maSlowLen, 'omitnan');
        yVisible = [TT.Close(idx); maFast(idx); maMid(idx); maSlow(idx)];
        yTopMin = min(yVisible, [], 'omitnan');
        yTopMax = max(yVisible, [], 'omitnan');
    otherwise
        yTopMin = min(TT.Low(idx), [], 'omitnan');
        yTopMax = max(TT.High(idx), [], 'omitnan');
end

if isfinite(yTopMin) && isfinite(yTopMax)
    d = yTopMax - yTopMin;
    if d == 0
        d = max(abs(yTopMax),1) * 0.01;
    end
    ylim(ax1, [yTopMin - 0.03*d, yTopMax + 0.03*d]);
end

yVolMax = max(TT.Volume(idx), [], 'omitnan');
if isfinite(yVolMax)
    d = max(yVolMax,1);
    ylim(ax2,[0, yVolMax + 0.05*d]);
end

if ~isempty(bottomFixedLims)
    ylim(ax3,bottomFixedLims);
else
    switch lower(bottomPanelMode)
        case 'macd'
            yBot = [macdLine(idx); macdSignal(idx); macdHist(idx)];
            yMin = min(yBot, [], 'omitnan');
            yMax = max(yBot, [], 'omitnan');
        case 'roc'
            yMin = min(rocVals(idx), [], 'omitnan');
            yMax = max(rocVals(idx), [], 'omitnan');
        otherwise
            return
    end

    if isfinite(yMin) && isfinite(yMax)
        d = yMax - yMin;
        if d == 0
            d = 1;
        end
        ylim(ax3,[yMin - 0.08*d, yMax + 0.08*d]);
    end
end
end

function addBottomIndicatorRows(hObj, bottomPanelMode, ...
    rsiVals, stochK, stochD, macdLine, macdSignal, macdHist, williamsR, rocVals)

switch lower(bottomPanelMode)
    case 'rsi'
        hObj.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('RSI', round(rsiVals,2));
    case 'stoch'
        hObj.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('%K', round(stochK,2));
        hObj.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('%D', round(stochD,2));
    case 'macd'
        hObj.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('MACD', round(macdLine,4));
        hObj.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Signal', round(macdSignal,4));
        hObj.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Hist', round(macdHist,4));
    case 'williamsr'
        hObj.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('Williams %R', round(williamsR,2));
    case 'roc'
        hObj.DataTipTemplate.DataTipRows(end+1) = dataTipTextRow('ROC (%)', round(rocVals,2));
    otherwise
        error('FSDA:getYahoo:WngInp','Unknown bottomPanelMode.');
end
end

function s = fmtVolume(v)
if v >= 1e9
    s = sprintf('%.2fB', v/1e9);
elseif v >= 1e6
    s = sprintf('%.2fM', v/1e6);
elseif v >= 1e3
    s = sprintf('%.1fK', v/1e3);
else
    s = sprintf('%.0f', v);
end
end

function showPanelExplanationCommand(topPanelMode, bottomPanelMode, ...
    maFastLen, maMidLen, maSlowLen, ...
    rsiLen, stochLen, stochSmooth, macdFastLen, macdSlowLen, macdSigLen, rocLen, wrLen, ...
    topRSI, lowRSI, topStoch, lowStoch, topWilliams, lowWilliams)

fprintf('\n');
fprintf('============================================================\n');
fprintf('EXPLANATION OF THE SELECTED PANELS\n');
fprintf('============================================================\n');

fprintf('\nTOP PANEL\n');
fprintf('------------------------------------------------------------\n');

switch lower(topPanelMode)
    case 'candle'
        fprintf('Mode: Candlestick chart\n');
        fprintf('  - Open, High, Low, Close for each bar.\n');
    case 'line'
        fprintf('Mode: Close price line\n');
        fprintf('  - Only the Close price is shown.\n');
    case 'ma'
        fprintf('Mode: Close + moving averages\n');
        fprintf('  - Fast MA  = %d\n', maFastLen);
        fprintf('  - Medium MA= %d\n', maMidLen);
        fprintf('  - Slow MA  = %d\n', maSlowLen);
end

fprintf('\nBOTTOM PANEL\n');
fprintf('------------------------------------------------------------\n');

switch lower(bottomPanelMode)
    case 'rsi'
        fprintf('Indicator: RSI\n');
        fprintf('  - WindowSize = %d\n', rsiLen);
        fprintf('  - Overbought threshold = %.0f\n', topRSI);
        fprintf('  - Oversold threshold = %.0f\n', lowRSI);
    case 'stoch'
        fprintf('Indicator: Stochastic Oscillator\n');
        fprintf('  - Lookback length = %d\n', stochLen);
        fprintf('  - Smoothing for %%D = %d\n', stochSmooth);
        fprintf('  - Upper threshold = %.0f\n', topStoch);
        fprintf('  - Lower threshold = %.0f\n', lowStoch);
    case 'macd'
        fprintf('Indicator: MACD\n');
        fprintf('  - Fast EMA   = %d\n', macdFastLen);
        fprintf('  - Slow EMA   = %d\n', macdSlowLen);
        fprintf('  - Signal EMA = %d\n', macdSigLen);
    case 'williamsr'
        fprintf('Indicator: Williams %%R\n');
        fprintf('  - Lookback length = %d\n', wrLen);
        fprintf('  - Upper threshold = %.0f\n', topWilliams);
        fprintf('  - Lower threshold = %.0f\n', lowWilliams);
    case 'roc'
        fprintf('Indicator: ROC\n');
        fprintf('  - Lag = %d\n', rocLen);
end

fprintf('\n============================================================\n\n');
end

function [intervalOut, isValid] = validateYahooRangeInterval(rangeIn, intervalIn, autoFix)
rangeIn = char(string(rangeIn));
validRanges = {'1d' '5d' '1mo' '3mo' '6mo' '1y' '2y' '5y' '10y' 'ytd' 'max'};

if ~ismember(rangeIn, validRanges)
    error('FSDA:getYahoo:WngInp', 'Unsupported range "%s".', rangeIn);
end

intervalIn = char(string(intervalIn));
validIntervals = {'1m','2m','5m','15m','30m','60m','90m','1h','1d','5d','1wk','1mo','3mo'};

if ~ismember(intervalIn, validIntervals)
    error('FSDA:getYahoo:WngInp', 'Unsupported interval "%s".', intervalIn);
end

switch rangeIn
    case '1d'
        allowed = {'1m','2m','5m','15m','30m','60m','90m','1h','1d'};
    case '5d'
        allowed = {'1m','2m','5m','15m','30m','60m','90m','1h','1d'};
    case '1mo'
        allowed = {'2m','5m','15m','30m','60m','90m','1h','1d','5d'};
    case '3mo'
        allowed = {'2m','5m','15m','30m','60m','90m','1h','1d','5d','1wk'};
    otherwise
        allowed = {'1d','5d','1wk','1mo','3mo'};
end

isValid = ismember(intervalIn, allowed);

if isValid
    intervalOut = intervalIn;
    return
end

if ~autoFix
    intervalOut = intervalIn;
    return
end

priority = {'1m','2m','5m','15m','30m','60m','90m','1h','1d','5d','1wk','1mo','3mo'};
idxAllowed = find(ismember(priority, allowed), 1, 'first');
intervalOut = priority{idxAllowed};
end

%FScategory:UTI-FIN
