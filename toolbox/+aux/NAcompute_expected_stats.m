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

zerop1=zeros(p,1);
zeropp=zeros(p);
T1 = zerop1;
T2 = zeropp;

for i = 1:n
    x_i = X(i, :);
    w_i = w(i);
    isnanxi=isnan(x_i);
    obs_ind = find(~isnanxi);     % observed indices
    mis_ind = find(isnanxi);      % missing indices

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

end