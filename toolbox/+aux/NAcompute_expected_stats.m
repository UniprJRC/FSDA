function [T1, T2] = NAcompute_expected_stats(X, mu, Sigma, w)
% NAcompute_expected_stats  E-step: compute expected sums T1 and T2 with missing data.
%
%
% Required input arguments:
%
%   X     - n x p data matrix with NaN for missing entries
%   mu    - p x 1 mean vector
%   Sigma - p x p covariance matrix
%
%  Optional input arguments:
%
%   w     - n x 1 weight vector (0/1 or nonnegative weights)
%
%
%  Output:
%
%   T1    - p x 1 weighted sum of expected x: sum_i w_i * E[x_i]
%   T2    - p x p weighted sum of expected outer products:
%           sum_i w_i * E[x_i x_i']
%
% Copyright 2008-2026.
% Written by FSDA team
%

%% Beginning of code

[n, p] = size(X);
if nargin < 4
    w = ones(n,1);
end
if p>5
    zerop1=zeros(p,1);
    zeropp=zeros(p);
    T1 = zerop1;
    T2 = zeropp;
    seqp=1:p;
    for i = 1:n
        x_i = X(i, :);
        w_i = w(i);
        isnanxi=isnan(x_i);
        % obs_ind = find(~isnanxi);     % observed indices
        % mis_ind = find(isnanxi);      % missing indices
        obs_ind = seqp(~isnanxi);     % observed indices
        mis_ind = seqp(isnanxi);      % missing indices

        if isempty(mis_ind)
            % Case 1: no missing values in the row
            x_exp = x_i(:);              % expected = observed
            C_exp = x_exp * x_exp';      % outer product x x'
        elseif isempty(obs_ind)
            % Case 2: all values missing
            x_exp = mu(:);               % expectation is the mean
            C_exp = Sigma + mu(:) * mu(:)'; % E[xx'] = Var + mu mu'
        else
            % Case 3: some missing and some observed
            % Partition mean and covariance
            mu_obs = mu(obs_ind);
            mu_mis = mu(mis_ind);

            Sigma_obs_obs = Sigma(obs_ind, obs_ind);
            Sigma_mis_obs = Sigma(mis_ind, obs_ind);
            Sigma_obs_mis = Sigma(obs_ind, mis_ind);
            Sigma_mis_mis = Sigma(mis_ind, mis_ind);

            % Regression coefficients for missing on observed
            % B = Sigma_mis_obs * inv(Sigma_obs_obs)
            B = Sigma_mis_obs / Sigma_obs_obs; % p_mis x p_obs

            % conditional expectation of missing given observed
            x_obs = x_i(obs_ind)';
            x_mis_exp = mu_mis + B * (x_obs - mu_obs);

            % conditional covariance of missing given observed
            Sigma_mis_cond = Sigma_mis_mis - B * Sigma_obs_mis;

            % full expected vector
            x_exp = zerop1;
            x_exp(obs_ind) = x_obs;
            x_exp(mis_ind) = x_mis_exp;

            % expected outer product matrix E[x x']
            C_exp = zeropp;
            % observed-observed block = observed outer product
            C_exp(obs_ind, obs_ind) = x_obs * x_obs';
            % observed-missing and missing-observed cross-products (use expectation for missing)
            C_exp(obs_ind, mis_ind) = x_obs * x_mis_exp';
            C_exp(mis_ind, obs_ind) = (C_exp(obs_ind, mis_ind))';
            % missing-missing block = conditional covariance + mean*mean'
            C_exp(mis_ind, mis_ind) = Sigma_mis_cond + x_mis_exp * x_mis_exp';
        end

        % accumulate weighted sums
        T1 = T1 + w_i * x_exp;
        T2 = T2 + w_i * C_exp;
    end
else
    mu = mu(:);
    T1 = zeros(p,1, 'like', X);
    T2 = zeros(p,p, 'like', X);

    % Missingness pattern
    M = isnan(X);

    % Group equal NaN patterns
    [patterns, ~, G] = unique(M, 'rows', 'stable');
    nPat = size(patterns,1);

    seqp = 1:p;
    muMuT = mu * mu.';

    for g = 1:nPat
        idx = (G == g);
        wg  = w(idx);
        sw  = sum(wg);

        if sw == 0
            continue
        end

        pat = patterns(g,:);
        obs_ind = seqp(~pat);
        mis_ind = seqp(pat);

        % Case 1: no missing values
        if isempty(mis_ind)
            Xg = X(idx,:);                       % ng x p
            T1 = T1 + Xg.' * wg;                % p x 1
            T2 = T2 + Xg.' * (Xg .* wg);        % p x p
            continue
        end

        % Case 2: all values missing
        if isempty(obs_ind)
            T1 = T1 + sw * mu;
            T2 = T2 + sw * (Sigma + muMuT);
            continue
        end

        % Case 3: partial missingness
        mu_obs = mu(obs_ind);
        mu_mis = mu(mis_ind);

        Sigma_obs_obs = Sigma(obs_ind, obs_ind);
        Sigma_mis_obs = Sigma(mis_ind, obs_ind);
        Sigma_obs_mis = Sigma(obs_ind, mis_ind);
        Sigma_mis_mis = Sigma(mis_ind, mis_ind);

        % Pattern-dependent quantities: computed once for the whole group
        B = Sigma_mis_obs / Sigma_obs_obs;              % nMis x nObs
        Sigma_mis_cond = Sigma_mis_mis - B * Sigma_obs_mis;

        Xobs = X(idx, obs_ind);                         % ng x nObs

        % E[x_mis | x_obs] for all rows in the group
        Xmis_exp = (mu_mis.' + (Xobs - mu_obs.') * B.');   % ng x nMis

        % ---- T1 ----
        T1(obs_ind) = T1(obs_ind) + Xobs.' * wg;
        T1(mis_ind) = T1(mis_ind) + Xmis_exp.' * wg;

        % ---- T2 ----
        % observed-observed
        T2(obs_ind, obs_ind) = T2(obs_ind, obs_ind) + Xobs.' * (Xobs .* wg);

        % observed-missing
        OM = Xobs.' * (Xmis_exp .* wg);
        T2(obs_ind, mis_ind) = T2(obs_ind, mis_ind) + OM;
        T2(mis_ind, obs_ind) = T2(mis_ind, obs_ind) + OM.';

        % missing-missing
        T2(mis_ind, mis_ind) = T2(mis_ind, mis_ind) ...
            + Xmis_exp.' * (Xmis_exp .* wg) ...
            + sw * Sigma_mis_cond;
    end
end

end