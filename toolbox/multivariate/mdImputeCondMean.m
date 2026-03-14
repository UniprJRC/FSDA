function Yimp = mdImputeCondMean(Y, mu, Sigma)
% Replace NaNs in Y by E[y_mis | y_obs, mu, Sigma] under MVN(mu,Sigma)
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
% Copyright  2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('mdImputeCondMean')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit


% Examples:

%{

    % Example of imputation using conditional mean.
    p=10; n=20;
    Yfull=randn(n,p);
    MU=median(Yfull); 
    Sigma=eye(p);
    Y=Yfull;
    Y(randsample(n*p,10,false))=NaN;
    Yinp=mdImputeCondMean(Y,MU,Sigma);
%}

%% Beginning of code

n = size(Y,1);
Yimp = Y;
% Column vector.
mu=mu(:);

for i = 1:n
    oi = ~isnan(Y(i,:));   % observed mask (1 x p)
    mi = ~oi;              % missing mask

    if ~any(mi)
        continue; % complete row
    end

    obsIdx = find(oi);
    misIdx = find(mi);

    y_obs = Y(i, obsIdx)';         % (#obs x 1)
    mu_obs = mu(obsIdx);           % (#obs x 1)
    mu_mis = mu(misIdx);           % (#mis x 1)

    S_oo = Sigma(obsIdx, obsIdx);
    S_mo = Sigma(misIdx, obsIdx);

    % Numerically stable solve: S_oo \ (y_obs - mu_obs)
    delta = S_oo \ (y_obs - mu_obs);

    y_mis_hat = mu_mis + S_mo * delta; % (#mis x 1)

    Yimp(i, misIdx) = y_mis_hat';
end

end