function d2_adj = mdPartialMD2full(d2, p, pobs, varargin)
% mdPartialMD2full  Rescale partial squared Mahalanobis distances to the full-dimensional scale.
%
%<a href="matlab: docsearchFS('mdPartialMD2full')">Link to the help function</a>
%
% Required input arguments:
%
%   d2    : squared partial distances. Vector.
%           Vector of length n containing the squared Mahalanobis distances
%           computed using only the observed variables for each row.
%           Data Types - single | double
%
%   p     : full dimension. Positive integer scalar.
%           p is the number of variables in the original data matrix.
%           Data Types - single | double
%
%   pobs  : number of observed variables. Vector.
%           Vector of length n containing the number of observed variables
%           for each row.
%           Data Types - single | double
%
% Optional input arguments:
%
%   method : method used to rescale the distances. String scalar or char vector.
%            Possible values are.
%
%            'pri'      : principled EM rescaling (default), 
%                         d2_partial + (p - pobs).
%
%            'expScale' : expectation scaling, 
%                         d2_partial * (p / pobs).
%
%            'zMap'     : standardization mapping,
%                         p + sqrt(2*p) * ((d2_partial - pobs) ./ sqrt(2*pobs)).
%
%            'detMap'   : determinant-based rescaling,
%                         d2_partial * (p / pobs) * (g_full / g_obs)
%                         This method requires input option 'Y' and input
%                         option 'Sigma'.
%
%            'chiMap'   : chi-square quantile mapping. Use the cdf and
%                         inverse of the cdf of Chi2 distribution.
%
%            'betaMap'  : Beta quantile mapping. Use the cdf and
%                         inverse of the cdf of Beta distribution.
%            Example - 'method','chiMap'
%            Data Types - string scalar | char vector
%
%
%   Y      : original data matrix with possible missing values. Matrix.
%            n x p data matrix. This input is required only if
%            method='detMap'.
%            Example - 'Y',Y
%            Data Types - single | double
%
%   Sigma  : covariance matrix in the full space. Matrix.
%            p x p covariance matrix. This input is required only if
%            method='detMap'.
%            Example - 'Sigma',Sigma
%            Data Types - single | double
%
% Output:
%
%   d2_adj : adjusted squared distances. Vector.
%            Column vector of length n containing the adjusted squared
%            Mahalanobis distances on the full p-dimensional scale.
%            If pobs(i)=0 or d2(i) is NaN, the corresponding output is NaN.
%
%
% More About:
%
%   The available methods are:
%
%   1) 'impMD'    : MD on EM-imputed data
%   2) 'pri'      : principled EM rescaling
%   3) 'expScale' : expectation scaling
%   4) 'zMap'     : standardization mapping
%   5) 'detMap'   : determinant-based rescaling
%   6) 'chiMap'   : chi-square mapping
%   7) 'betaMap'  : Beta mapping
%
%
% References:
%
%   Little, R. J. A., & Rubin, D. B. (2020). Statistical Analysis with
%   Missing Data (3rd ed.). Hoboken, NJ: John Wiley & Sons.
%
%   Wilks, S. S. (1962). Mathematical Statistics. John Wiley & Sons,
%   New York.
%
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdPartialMD2full')">Link to the help page for this function</a>
% 
%$LastChangedDate::                      $: Date of the last commit
%
%
% See also: mdPartialMD, mdEM, mdImputeCondMean
%
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Example 1: default principal mapping.
    p = 5;
    n = 100;
    rho = 0.9;
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);
    missRate = 0.01;
    Yfull = randn(n,p) * R;
    missMask = rand(n,p) < missRate;
    Y = Yfull;
    Y(missMask) = NaN;
    mu = mean(Yfull)';
    [d2partial,pobs] = mdPartialMD(Y,mu,Sigma);
    d2_adj = mdPartialMD2full(d2partial,p,pobs);
%}
%
%{
    % Example 2: chi-square mapping.
    p = 5;
    n = 100;
    rho = 0.9;
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);
    missRate = 0.01;
    Yfull = randn(n,p) * R;
    missMask = rand(n,p) < missRate;
    Y = Yfull;
    Y(missMask) = NaN;
    mu = mean(Yfull)';
    [d2partial,pobs] = mdPartialMD(Y,mu,Sigma);
    d2_adj = mdPartialMD2full(d2partial,p,pobs,'method','chiMap');
%}
%
%{
    % Example 3: principled EM rescaling.
    p = 5;
    n = 100;
    rho = 0.9;
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);
    missRate = 0.01;
    Yfull = randn(n,p) * R;
    missMask = rand(n,p) < missRate;
    Y = Yfull;
    Y(missMask) = NaN;
    mu = mean(Yfull)';
    [d2partial,pobs] = mdPartialMD(Y,mu,Sigma);
    d2_adj = mdPartialMD2full(d2partial,p,pobs,'method','pri');
%}
%
%{
    % Example 4: determinant-based rescaling.
    p = 5;
    n = 100;
    rho = 0.9;
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);
    missRate = 0.01;
    Yfull = randn(n,p) * R;
    missMask = rand(n,p) < missRate;
    Y = Yfull;
    Y(missMask) = NaN;
    mu = mean(Yfull)';
    outEM = mdEM(Y,'condmeanimp',true);
    SigHat = outEM.cov;
    [d2partial,pobs] = mdPartialMD(Y,mu,SigHat);
    d2_adj = mdPartialMD2full(d2partial,p,pobs,'method','detMap', ...
        'Y',Y,'Sigma',Sigma);
%}

%% Beginning of code

d2 = d2(:);
pobs = pobs(:);

if numel(d2) ~= numel(pobs)
    error('FSDA:mdPartialMD2full:WrongInputOpt', ...
        'd2 and pobs must have the same length');
end
method='pri';

options = struct('method',method,'Y',[],'Sigma',[]);

if nargin>3
    [varargin{:}] = convertStringsToChars(varargin{:});
    UserOptions = varargin(1:2:length(varargin));

    if ~isempty(UserOptions)
        if length(varargin) ~= 2*length(UserOptions)
            error('FSDA:mdPartialMD2full:WrongInputOpt', ...
                ['Number of supplied options is invalid. Probably values ', ...
                 'for some parameters are missing.']);
        end
        aux.chkoptions(options,UserOptions)
    end

    for i=1:2:length(varargin)
        options.(varargin{i}) = varargin{i+1};
    end
end

method = options.method;
if isstring(method)
    method = char(method);
end

if ~ischar(method)
    error('FSDA:mdPartialMD2full:WrongInputOpt', ...
        'Option ''method'' must be a string scalar or a char vector');
end

validMethods = {'pri','expScale','zMap','detMap','chiMap','betaMap'};
if ~any(strcmp(method,validMethods))
    error('FSDA:mdPartialMD2full:WrongInputOpt', ...
        ['Option ''method'' must be one of: ''pri'', ', ...
         '''expScale'', ''zMap'', ''detMap'', ''chiMap'', ''betaMap''.']);
end

n = numel(d2);
d2_adj = nan(n,1);

switch method

    case 'pri'
        idx = pobs>0 & ~isnan(d2);
        d2_adj(idx) = d2(idx) + (p - pobs(idx));

    case 'expScale'
        idx = pobs>0 & ~isnan(d2);
        d2_adj(idx) = d2(idx) .* (p ./ pobs(idx));

    case 'zMap'
        idx = pobs>0 & ~isnan(d2);
        z = (d2(idx) - pobs(idx)) ./ sqrt(2*pobs(idx));
        d2_adj(idx) = p + sqrt(2*p) .* z;

    case 'detMap'
        Y = options.Y;
        Sigma = options.Sigma;

        if isempty(Y) || isempty(Sigma)
            error('FSDA:mdPartialMD2full:WrongInputOpt', ...
                'Method ''detMap'' requires input options ''Y'' and ''Sigma''.');
        end

        if size(Y,1) ~= n
            error('FSDA:mdPartialMD2full:WrongInputOpt', ...
                'Input Y must have the same number of rows as length(d2).');
        end

        if size(Y,2) ~= p
            error('FSDA:mdPartialMD2full:WrongInputOpt', ...
                'Input Y must have p columns.');
        end

        if ~isequal(size(Sigma),[p p])
            error('FSDA:mdPartialMD2full:WrongInputOpt', ...
                'Input Sigma must be a p-by-p matrix.');
        end

        d2_adj = localDetMapFromPartial(Y,d2,pobs,Sigma,p);

    case 'chiMap'
        epsv = eps;
        for i = 1:n
            p1 = pobs(i);

            if isnan(d2(i)) || p1<=0
                d2_adj(i) = NaN;
                continue
            end

            c = chi2cdf(d2(i),p1);
            c = min(max(c,epsv),1-epsv);
            d2_adj(i) = chi2inv(c,p);
        end

    case 'betaMap'
        c1 = (n-1)^2 / n;
        epsv = eps;

        a2 = p/2;
        b2 = (n-p-1)/2;

        if b2 <= 0
            error('FSDA:mdPartialMD2full:WrongInputOpt', ...
                'For method ''betaMap'' it must be n > p + 1.');
        end

        for i=1:n
            p1 = pobs(i);

            if isnan(d2(i)) || p1<=0 || n<=p1+1
                d2_adj(i) = NaN;
                continue
            end

            a1 = p1/2;
            b1 = (n-p1-1)/2;

            u = d2(i)/c1;
            u = min(max(u,0),1-epsv);

            alpha = betacdf(u,a1,b1);
            alpha = min(max(alpha,epsv),1-epsv);

            v = betainv(alpha,a2,b2);
            d2_adj(i) = c1*v;
        end



end

end

%% Local functions

function d2_det = localDetMapFromPartial(Y, d2_part, poss, Sigma, p)
nloc = size(Y,1);
d2_det = nan(nloc,1);

logdet_full = localLogdetSPD(Sigma);
g_full = exp(logdet_full / p);

idx = poss>0 & ~isnan(d2_part);
ii = find(idx);

for j = 1:numel(ii)
    i = ii(j);
    obs = ~isnan(Y(i,:));
    pii = poss(i);

    Sobs = Sigma(obs,obs);
    logdet_obs = localLogdetSPD(Sobs);
    g_obs = exp(logdet_obs / pii);

    d2_det(i) = d2_part(i) * (p / pii) * (g_full / g_obs);
end

end

function val = localLogdetSPD(S)
R = chol(S);
val = 2*sum(log(diag(R)));
end

%FScategory:MULT-MissingData