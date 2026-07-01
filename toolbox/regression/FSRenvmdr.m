function [MDRenv] = FSRenvmdr(n,p,varargin)
%FSRenvmdr computes the theoretical envelopes of Minimum Deletion Residual outside subset during the search
%
% Optimized: replaces finv/tinv/norminv with direct betaincinv/erfcinv
% to eliminate Statistics Toolbox per-call overhead.
%
%<a href="matlab: docsearchFS('FSRenvmdr')">Link to the help function</a>
%
%  Required input arguments:
%
%    n : number of observations. Scalar.
%    p : number of explanatory variables (including the intercept if present). Scalar.
%
%  Optional input arguments:
%
%   init:       Search initialization. Scalar.
%   prob:       quantiles for which envelopes have to be computed. Vector.
%
%  Output:
%
%  MDRenv:      forward envelopes of mdr. Matrix.

% Copyright 2008-2025. Written by FSDA team
% Optimized by Sindhuja Parimalarangan, June 2026

%% Beginning of code

if ~isscalar(n) || isempty(n) || isnan(n)
    error('FSDA:FSRenvmdr:missingInputs','n must be scalar non empty and non missing!!');
end

if ~isscalar(p) || isempty(n) || isnan(p)
    error('FSDA:FSRenvmdr:missingInputs','p must be scalar non empty and non missing!!!');
end

if n<40
    inisearch=p+1;
else
    inisearch=min(3*p+1,floor(0.5*(n+p+1)));
end

if coder.target('MATLAB')
    prob=[0.01 0.5 0.99];
    options=struct('init',inisearch,'prob',prob);
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions=varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:FSRenvmdr:WrongInputOpt','Number of supplied options is invalid. Probably values for some parameters are missing.');
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

if m0>n-1
    error('FSDA:FSRenvmdr:TooLargen','Initial starting point of the search (m0=%d) is greater than n-1 (n-1=%d)', m0, n-1);
end

%% Envelopes generation

if size(prob,1)>1
    probf=prob';
else
    probf=prob;
end

lp = length(probf);
probf = 1-probf;

m =(m0:n-1)';
lm=length(m);
mm = repmat(m,1,lp);

% Direct finv via betaincinv (bypass Statistics Toolbox overhead)
d1 = 2*(n-mm);
d2 = 2*(mm+1);
x = betaincinv(repmat(probf,lm,1), d1/2, d2/2);
quant = (d2 ./ d1) .* x ./ (1 - x);

% from the equivalence with Incomplete beta distribution.
q=(mm+1)./(mm+1+(n-mm).*(quant));

% Direct tinv via betaincinv: tinv(p,nu) for p > 0.5
% tinv(0.5*(1+q), mm-p) = sqrt(nu * (1/betaincinv(2*(1-p), nu/2, 0.5) - 1))
% where p = 0.5*(1+q), so 2*(1-p) = 1-q
nu = mm - p;
xb = betaincinv(1-q, nu/2, 0.5);
MinSca = sqrt(nu .* (1./xb - 1));

% Direct norminv via erfcinv: norminv(p) = -sqrt(2)*erfcinv(2*p)
a = -sqrt(2) * erfcinv(2 * 0.5*(1+mm/n));
corr=1-2*(n./mm).*a.*exp(-0.5*a.^2)/sqrt(2*pi);

MDRenv=[m MinSca./sqrt(corr)];

end
%FScategory:REG-Regression
