function [d2, pobs] = mdPartialMD(Y, mu, Sigma)
% mdPartialMD computes squared Mahalanobis distances using only observed entries.
%
%
%<a href="matlab: docsearchFS('mdPartialMD')">Link to the help function</a>
%
% The function computes for each row i the (squared) Mahalanobis distance:
%   d2(i) = (x_obs - mu_obs)' * inv(Sigma_obs_obs) * (x_obs - mu_obs)
% where obs are the indices of non-missing entries in row i.
%
%
% Required input arguments:
%
%
%   Y      : data matrix. 2D array.
%            n x p data matrix (rows may contain NaN for missing)
%   mu     : location estimate. Vector. Vector of length p containing
%            estimate of location 
% Sigma    : scatter estimate. Matrix. p x p estimated  covariance matrix.
%
%  Optional input arguments:
%
%  Output:
%
%   d2     : squared Mahalanobis distances (computed on
%           observed dims). n x 1 vector. If a row has zero observed
%           entries, d2(i) is set to NaN.
%   pobs   : number of observed variables  for each row. Vector. 
%           n x 1 vector with number of observed variable for each row of
%           input data matrix X. If a row has zero observed entries
%           poss(i)=0.
%
%
% See also: mdEM, mdTEM.m, mdImputeCondMean.m, mdPartialMD2full
%
%
% References:
%
% Little, R. J. A., & Rubin, D. B. (2019). Statistical Analysis with
% Missing Data (3rd ed.). Hoboken, NJ: John Wiley & Sons.
% van Buuren, S. (2018). Flexible Imputation of Missing Data (2nd ed.).
% Boca Raton, FL: Chapman & Hall/CRC (Taylor & Francis Group).
% Templ, M. (2023). Visualization and Imputation of Missing Values: With
% Applications in R. Cham, Switzerland: Springer Nature. 
%
%
% Copyright 2008-2026.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('mdPartialMD')">Link to the help page for this function</a>
%

%{
    % Example with just one output argument.
    p = 5;                % number of variables
    n = 100;             % number of observations
    rho = 0.9;            % target pairwise correlation (0<rho<1)
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);      % upper-triangular such that Sigma = R'*R
    missRate = 0.01;     % MCAR missing probability per entry
    
    % Generate samples ~ N(0,Sigma)
    Yfull = randn(n,p) * R;   % Strong positive correlation between the vars
    missMask = rand(n,p) < missRate;
    mu=mean(Yfull)';
    S=cov(Yfull);
    Y=Yfull;
    Y(missMask) = NaN;
    d2fast=mdPartialMD(Y,mu,Sigma);
%}

%{
    % Example with two output arguments.
    p = 5;                % number of variables
    n = 100;             % number of observations
    rho = 0.9;            % target pairwise correlation (0<rho<1)
    Sigma = (1-rho)*eye(p) + rho*ones(p);
    R = chol(Sigma);      % upper-triangular such that Sigma = R'*R
    missRate = 0.01;     % MCAR missing probability per entry
    
    % Generate samples ~ N(0,Sigma)
    Yfull = randn(n,p) * R;   % Strong positive correlation between the vars
    missMask = rand(n,p) < missRate;
    mu=mean(Yfull)';
    S=cov(Yfull);
    Y=Yfull;
    Y(missMask) = NaN;
    [d2fast,pobs]=mdPartialMD(Y,mu,Sigma);
%}

%% Beginning of code

[n,p] = size(Y);

d2 = nan(n,1);
pobs = zeros(n,1);
mu=mu(:);

if n>=300 && p<=5
    isnanX=isnan(Y);
    patt=unique(isnanX,"rows");

    for i = 1:size(patt,1)
        % Select all rows of matrix X which have missing patter
        patti=patt(i,:);

        seli=all(isnanX==patti,2);
        pattiNot=~patti;
        Y_i = Y(seli,pattiNot);

        pobs(seli)=p-sum(patti);

        % Extract observed parts of x, mu, and Sigma
        mu_obs = mu(pattiNot);          % p1 x 1
        Sigma_obs = Sigma(pattiNot, pattiNot); % p1 x p1

        % squared Mahalanobis distance
        d2(seli)=mahalFS(Y_i,mu_obs',Sigma_obs);
    end

else

    for i = 1:n
        x_i = Y(i, :);                 % 1 x p
        obs_ind = find(~isnan(x_i));   % observed indices for row i
        p1 = numel(obs_ind);           % number observed

        pobs(i) = p1;

        if p1 == 0
            % No observed data => cannot compute Mahalanobis
            d2(i) = NaN;
            continue
        end

        % Extract observed parts of x, mu, and Sigma
        x_obs = x_i(obs_ind)';         % p1 x 1
        mu_obs = mu(obs_ind);          % p1 x 1
        Sigma_obs = Sigma(obs_ind, obs_ind); % p1 x p1

        % If p1 == 1, Sigma_obs is scalar; guard inversion accordingly
        % Use stable solver rather than explicit inv()
        diff = x_obs - mu_obs;
        % Solve Sigma_obs * y = diff  => y = Sigma_obs \ diff
        % Then d2 = diff' * (Sigma_obs \ diff)
        y = Sigma_obs \ diff;
        d2(i) = diff' * y;             % squared Mahalanobis distance
    end
end