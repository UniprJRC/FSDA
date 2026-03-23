function Yimp = mdImputeStochastic(Y, mu, Sigma)
% Replace NaNs in Y by random draws from
% Y_mis | Y_obs ~ N(mu_cond, Sigma_cond) under MVN(mu,Sigma)
%
%<a href="matlab: docsearchFS('mdImputeCondMean')">Link to the help function</a>
%
%  Required input arguments:
%
% Y :           Input data. Matrix. 
%               n x v data matrix; n observations and v variables. Rows of
%               Y represent observations, and columns represent variables.
%               Y contains a seres of missing values (NaN's) 
%                Data Types - single|double
%        mu :   Centroid. Vector.  vector of length v, containing centroid to use
%      Sigma:   Covariance matrix. Matrix. v x v matrix containing covariance matrix which must be used
%       
%
%  Optional input arguments:
%
%  Output:
%
%    Yimp :   matrix with imputed values. Matrix of size nxv. Matrix Yinp
%            does not contain missing values.
%
% See also: mdEM
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
% Notes:
% - This performs single stochastic imputation.
% - Re-running the function gives different imputations unless rng is fixed.
% - For proper multiple imputation, repeat this several times and analyze
%   the completed datasets separately.
%
% Copyright  2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mdImputeCondMean')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


% Examples:

%{

    % Example of stochastic imputation, conditional mu and sigma.
    p=10; n=20;
    Yfull=randn(n,p);
    MU=median(Yfull); 
    Sigma=eye(p);
    Y=Yfull;
    Y(randsample(n*p,10,false))=NaN;
    Yinp=mdImputeStochastic(Y,MU,Sigma);
%}

%% Beginning of code



n = size(Y,1);
Yimp = Y;
mu = mu(:);   % force column vector

for i = 1:n
    oi = ~isnan(Y(i,:));   % observed entries
    mi = ~oi;              % missing entries

    if ~any(mi)
        continue;   % row already complete
    end

    obsIdx = find(oi);
    misIdx = find(mi);

    y_obs  = Y(i, obsIdx)';          % (#obs x 1)
    mu_obs = mu(obsIdx);             % (#obs x 1)
    mu_mis = mu(misIdx);             % (#mis x 1)

    S_oo = Sigma(obsIdx, obsIdx);    % observed-observed block
    S_mo = Sigma(misIdx, obsIdx);    % missing-observed block
    S_om = Sigma(obsIdx, misIdx);    % observed-missing block
    S_mm = Sigma(misIdx, misIdx);    % missing-missing block

    % Conditional mean:
    % mu_cond = mu_mis + S_mo * inv(S_oo) * (y_obs - mu_obs)
    delta   = S_oo \ (y_obs - mu_obs);
    mu_cond = mu_mis + S_mo * delta;

    % Conditional covariance:
    % Sigma_cond = S_mm - S_mo * inv(S_oo) * S_om
    Sigma_cond = S_mm - S_mo * (S_oo \ S_om);

    % Force symmetry (helps against numerical roundoff)
    Sigma_cond = (Sigma_cond + Sigma_cond') / 2;

    % Draw one random vector from N(mu_cond, Sigma_cond)
    [L,p] = chol(Sigma_cond, 'lower');

    if p == 0
        % Positive definite case
        y_mis_draw = mu_cond + L * randn(length(misIdx),1);
    else
        % Fallback: eigenvalue cleanup for semidefinite / near-singular case
        [V,D] = eig(Sigma_cond);
        d = diag(D);
        d(d < 0) = 0;   % remove tiny negative eigenvalues from roundoff
        A = V * diag(sqrt(d));
        y_mis_draw = mu_cond + A * randn(length(misIdx),1);
    end

    % Put imputed values back into the row
    Yimp(i, misIdx) = y_mis_draw';
end

end