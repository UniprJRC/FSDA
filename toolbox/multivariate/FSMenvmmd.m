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
% n :           Number of observations. Scalar.
% v :           Number of variables. Scalar.
%
% Optional input arguments:
%
% init :       Point where to start monitoring. Scalar.
% prob:        quantiles for which envelopes have to be computed. Vector.
% scaled:     If true, no consistency factor applied. Boolean.
%
% Output:
%
%  MMDenv=      n-m0+1 x length(prob)+1 columns containing the envelopes.

% Copyright 2008-2025. Written by FSDA team
% Optimized by Sindhuja Parimalarangan, June 2026

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
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSMenvmmd:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
        end
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

% Compute Tallis correction factor
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
