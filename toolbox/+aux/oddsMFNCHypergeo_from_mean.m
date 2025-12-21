function odds = oddsMFNCHypergeo_from_mean(mu, m, n)
% Estimate odds from mean from multivariate Fisher non central hypergeometric
% distribution
% MATLAB translation of the C++ oddsMFNCHypergeo_from_mean routine.
colors = numel(mu);
if colors < 2
    error("FSDA:oddsMFNCHypergeo_from_mean","Number of colors too small");
end
if length(mu)~=colors
    error("FSDA:oddsMFNCHypergeo_from_mean","Length of vectors mu and m must be the same");
end

% Initialize output
odds = nan(colors,1);

% Just in case mu and m are row vectors

mu = mu(:);
m  = m(:);
colors = numel(m);


%% Validate  parameters
if n < 0 || any(m < 0)
    error("FSDA:oddsMFNCHypergeo_from_mean","n and m must not be negative")
end

if length(mu)~=length(m)
    % warning flag kept for compatibility
    warning("FSDA:oddsMFNCHypergeo_from_mean","Length of vectors mu and m must be the same");
end


%% Compute N = sum(m) and check sum(m)
N = sum(m);  % N = total number of balls in the urn
sum_mu=sum(mu);

if abs(sum_mu-n)/n>0.1
    error("FSDA:oddsMFNCHypergeo_from_mean","sum of means must be equal to n")
end
if n > N
    error("FSDA:oddsMFNCHypergeo_from_mean","n > sum(m): Taking more items than there are")
end

% Choose reference: max distance from bounds
c0 = 1; xd0 = 0;
for j = 1:colors
    x1 = m(j) + n - N; if x1 < 0, x1 = 0; end
    x2 = n; if x2 > m(j), x2 = m(j); end
    xd1 = mu(j) - x1;
    xd2 = x2 - mu(j);
    xd1 = min(xd1, xd2);
    if xd1 > xd0
        xd0 = xd1;
        c0 = j;
    end
end

%% If all odds undetermined
if xd0 == 0
    error("FSDA:oddsMFNCHypergeo_from_mean","all odds are indetermined")
end

odds(c0) = 1;

for j = 1:colors
    if j == c0, continue; end

    x1 = m(j) + n - N; if x1 < 0, x1 = 0; end
    x2 = n; if x2 > m(j), x2 = m(j); end
    muj = mu(j);

    if x1 == x2
        odds(j) = NaN;
        warning("FSDA:oddsMFNCHypergeo_from_mean","missing values in odds")
        continue
    end

    if muj <= x1
        if muj == x1
            odds(j) = 0.0;
            warning("FSDA:oddsMFNCHypergeo_from_mean","0 values in  odds")
        else
            odds(j) = NaN;
            warning("FSDA:oddsMFNCHypergeo_from_mean","NaN values in  odds")
        end
        continue
    end

    if muj >= x2
        if muj == x2
            odds(j) = Inf;
            warning("FSDA:oddsMFNCHypergeo_from_mean","Inf values in  odds")
        else
            odds(j) = NaN;
            warning("FSDA:oddsMFNCHypergeo_from_mean","NaN values in  odds")
        end
        continue
    end


    odds(j) = muj * (m(c0) - mu(c0)) / (mu(c0) * (m(j) - muj));
end
end