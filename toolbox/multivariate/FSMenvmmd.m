function [MMDenv] = FSMenvmmd(n,v,varargin)
%FSMenvmmd computes the theoretical envelopes of Minimum MD outside subset during the search
%
% Optimized: replaces finv/chi2inv/chi2cdf with betaincinv/gammaincinv/gammainc
% to eliminate Statistics Toolbox per-call overhead.
%
%<a href="matlab: docsearchFS('FSMenvmmd')">Link to the help function</a>
%
% Required input arguments:
%
% n :           Number of observations. Scalar. Number of observations.
%               Data Types - single | double
% v :           Number of variables. Scalar. Number of variables.
%               Data Types - single | double
%
% Optional input arguments:
%
% init :       Point where to start monitoring required diagnostics. Scalar.
%              Note that if bsb is supplied, init>=length(bsb). If init is not
%              specified it will be set equal to floor(n*0.6).
%                 Example - 'init',50
%                 Data Types - double
%
% prob:        quantiles for which envelopes have
%               to be computed. Vector. Vector containing 1 x k elements .
%               The default is to produce 1 per cent, 50 per cent and 99 per cent envelopes.
%                 Example - 'prob',[0.05 0.95]
%                 Data Types - double
%
%   scaled:  It indicates how to compute the envelopes. Boolean.
%               If scaled=true0 the envelopes are produced for
%               scaled Mahalanobis distances (no consistency factor is
%               applied) else the traditional consistency factor is applied
%               (this is the default)
%                 Example - 'scaled',false
%                 Data Types - logical
%
%
% Output:
%
%  MMDenv=      n-m0+1 x length(prob)+1 columns containing the envelopes
%               for the requested quantiles.
%               1st col = fwd search index from m0 to n-1;
%               2nd col = envelope for quantile prob[1];
%               3rd col = envelope for quantile prob[2];
%               ...;
%               (k+1) col = envelope for quantile prob[k].
%
% See also FSMenvmmd.m, FSM.m
%
% References:
%
% Riani, M., Atkinson, A.C. and Cerioli, A. (2009), Finding an unknown
% number of multivariate outliers, "Journal of the Royal Statistical
% Society Series B", Vol. 71, pp. 201-221.
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('FSMenvmmd')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%
%{
    %% FSMenvmmd with all default options.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent envelopes based on 10000
    % observations and 5 explanatory variables.
    MMDenv=FSMenvmmd(10000,5);
    plot(MMDenv(:,1),MMDenv(:,2:end))
%}

%{
    %% FSMenvmmd with optional arguments.
    % Example of creation of 1 per cent, 50 per cent and 99 per cent envelopes based on 10000
    % observations and 5 explanatory variables. The envelopes are produced for
    % scaled Mahalanobis distances (no consistency factor is applied)
    MMDenv=FSMenvmmd(10000,5,'scaled',1);
    plot(MMDenv(:,1),MMDenv(:,2:end))
%}

%{
    %% Order statistics and simulations envelopes .
    % In this example we compare the accuracy of the envelopes computed with
    % order statistics with those which come from simulations.

    % Fix a seed
    state=1000;

    mtstream = RandStream('shr3cong','Seed',state);
    RandStream.setGlobalStream(mtstream);
    defaultStream = RandStream.getGlobalStream();
    reset(defaultStream)

    % If you run this example in a version older than 7.9 replace the previous
    % four lines with
    % randn('state', 1000);
    n=200;
    p=3;


    init=25;
    nsimul=1000;
    mmdStore=zeros(n-init,nsimul);

    for j=1:nsimul
        Y=randn(n,p);
        [fre]=unibiv(Y);
        %create an initial subset with the 20 observations with the lowest
        %Mahalanobis Distance
        fre=sortrows(fre,4);
        bs=fre(1:25,1);
        mmd = FSMmmd(Y,bs,'init',init);
        mmdStore(:,j)=mmd(:,2);
    end

    % Sort rows of matrix mmdStore
    mmdStore=sort(mmdStore,2);

    % Create figure which compares empirical and theoretical forward envelopes
    % for minimum deletion residual
    figure;
    hold('on');
    quant=[0.01 0.5 0.99];
    sel=round(nsimul*quant);
    % Plot lines of empirical quantiles
    line(mmd(:,1),mmdStore(:,sel),'LineStyle','--','Color','g');
    % Plots lines of theoretical quantiles using order statistics
    mmdT=FSMenvmmd(n,p,'init',init);
    line(mmdT(:,1),mmdT(:,2:4),'LineStyle','-','Color','r');
    xlabel('Subset size m');
%}
%% Beginning of code

if ~isscalar(n) || isempty(n) || isnan(n)
    error('FSDA:FSMenvmmd:Wrongn','n must be scalar non empty and non missing!!');
end

if ~isscalar(v) || isempty(n) || isnan(v)
    error('FSDA:FSMenvmmd:Wrongv','v must be scalar non empty and non missing!!!');
end

inisearch=floor(n*0.6);

prob=[0.01 0.5 0.99];
scaled=false;
options=struct('init',inisearch,'prob',prob,'scaled',scaled);

if coder.target('MATLAB')
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        % Check if number of supplied options is valid
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSMenvmmd:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
        % Check if user options are valid options
        aux.chkoptions(options,UserOptions)
    end
end

if nargin>2
    for i=1:2:length(varargin)
        options.(varargin{i})=varargin{i+1};
    end
end

m0=options.init;
prob=options.prob;
scaled=options.scaled;

if m0>n-1
    if coder.target('MATLAB')
        error('FSDA:FSMenvmmd:WrongM0',['Initial starting point of the search (m0=' num2str(m0) ') is greater than n-1(n-1=' num2str(n-1) ')']);
    else
        error('FSDA:FSMenvmmd:WrongM0','Initial starting point of the search is greater than n-1');
    end
end

%% Envelopes generation

if size(prob,1)>1
    prob=prob';
end

lp = length(prob);
prob = 1-prob;

m =(m0:n-1)';
lm=length(m);
mm = repmat(m,1,lp);

% CODE BEFORE MATHWORKS OPTIMIZATION
% quant=finv(repmat(prob,lm,1),2*(n-mm),2*(mm+1));

% Direct finv via betaincinv (replaces finv from Statistics Toolbox)
d1 = 2*(n-mm);
d2 = 2*(mm+1);
x = betaincinv(repmat(prob,lm,1), d1/2, d2/2);
quant = (d2 ./ d1) .* x ./ (1 - x);

% from the equivalence with the incomplete beta distribution.
q=(mm+1)./(mm+1+(n-mm).*(quant));

cor=v*((mm+1)./mm).*(mm-1)./(mm-v);

% Direct finv for second call: finv(q, v, mm-v) via betaincinv
x2 = betaincinv(q, v/2, (mm-v)/2);
fval = ((mm-v)/v) .* x2 ./ (1 - x2);
MinSca = sqrt(cor .* fval);

% Compute Tallis correction factor based on the chi^2 distribution
% mm/n is the percentage of observations inside subset if scaled is not equal to 1.

if scaled == true
    corr = 1;
else
    % chi2inv(mm/n, v) via gammaincinv
    a = 2 * gammaincinv(mm/n, v/2);
    % chi2cdf(a, v+2) via gammainc
    corr = (n./mm) .* gammainc(a/2, (v+2)/2);
end

MMDenv=[m MinSca./sqrt(corr)];

end
%FScategory:MULT-Multivariate
