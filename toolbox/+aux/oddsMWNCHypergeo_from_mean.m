function odds = oddsMWNCHypergeo_from_mean(mu, m, n)
% Estimate odds from mean mu using Manly's approximation.
% MATLAB translation of the C++ oddsMWNCHypergeo_from_mean routine.

colors = numel(mu);
if colors < 2 
    error("FSDA:oddsMWNCHypergeo_from_mean","Number of colors too small");
end
if length(mu)~=colors
    error("FSDA:oddsMWNCHypergeo_from_mean","Length of vectors mu and m must be the same");
end

% Initialize output
odds = nan(colors,1);

% Just in case mu and m are row vectors
mu = mu(:);
m  = m(:);

%% Validate scalar parameters
if n < 0
    error("FSDA:oddsMWNCHypergeo_from_mean","n must be greater than 0")
end

if length(mu)~=length(m)
    % warning flag kept for compatibility
    warning("FSDA:oddsMWNCHypergeo_from_mean","Length of vectors mu and m must be the same");
end

%% Compute N = sum(m) and check m
N = sum(m);  % N = total number of balls in the urn
sum_mu=sum(mu);

if abs(sum_mu-n)/n>0.1
    error("FSDA:oddsMWNCHypergeo_from_mean","sum of means must be equal to n")
end
if n > N
    error("FSDA:oddsMWNCHypergeo_from_mean","n > sum(m): Taking more items than there are")
end


%% Find reference color
c0 = 1;
xd0 = 0.0;

for j = 1:colors
    x1 = m(j) + n - N;
    if x1 < 0, x1 = 0; end

    x2 = n;
    if x2 > m(j), x2 = m(j); end

    xd1 = mu(j) - x1;
    xd2 = x2 - mu(j);
    xd = min(xd1, xd2);

    if xd > xd0
        xd0 = xd;
        c0 = j;
    end
end

%% If all odds undetermined
if xd0 == 0
    error("FSDA:oddsMWNCHypergeo_from_mean","all odds are indetermined")
end

%% Set reference
odds(c0) = 1.0;

% Denominator
denom = log(1.0 - mu(c0) / m(c0));

%% Compute odds for other colors
for j = 1:colors
    if j == c0
        continue
    end

    x1 = m(j) + n - N;
    if x1 < 0, x1 = 0; end

    x2 = n;
    if x2 > m(j), x2 = m(j); end

    if x1 == x2
        odds(j) = NaN;
        warning("FSDA:oddsMWNCHypergeo_from_mean","missing values in odds")
        continue
    end

    muj = mu(j);

    if muj <= x1
        if muj == x1
            odds(j) = 0.0;
            warning("FSDA:oddsMWNCHypergeo_from_mean","0 values in  odds")
        else
            odds(j) = NaN;
            warning("FSDA:oddsMWNCHypergeo_from_mean","NaN values in  odds")
        end
        continue
    end

    if muj >= x2
        if muj == x2
            odds(j) = Inf;
            warning("FSDA:oddsMWNCHypergeo_from_mean","Inf values in  odds")
        else
            odds(j) = NaN;
            warning("FSDA:oddsMWNCHypergeo_from_mean","NaN values in  odds")
        end
        continue
    end

    % Manly odds ratio
    numer = log(1.0 - muj / m(j));
    odds(j) = numer / denom;
end

end

