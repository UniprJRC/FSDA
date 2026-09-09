function logf = logvMFpdfFS(Y, mu, kappa)
%logvMFpdfFS computes the log-density of the von Mises-Fisher distribution
%
%<a href="matlab: docsearchFS('logvMFpdfFS')">Link to the help function</a>
%
%   logvMFpdfFS computes, for an n-by-p matrix Y of unit vectors (rows
%   are points on the unit hypersphere S^(p-1) in R^p), the log-density of
%   a single von Mises-Fisher component with mean direction mu and
%   concentration kappa:
%
%       f(y; mu, kappa) = Cp(kappa) * exp(kappa * mu' * y)
%
%       Cp(kappa) = kappa^(p/2-1) / ( (2*pi)^(p/2) * I_(p/2-1)(kappa) )
%
%   where I_nu is the modified Bessel function of the first kind of order
%   nu. The normalizing constant is evaluated entirely in log-space using
%   the exponentially-scaled Bessel function (besseli(...,1)) to avoid
%   overflow for even moderate kappa or dimension p; this is the vMF
%   counterpart of the log-space evaluation performed by logmvnpdfFS for
%   the Gaussian density used in tclust.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-p matrix. Each row of Y is
%               assumed to be a unit vector (rows are NOT renormalized
%               inside this function: callers -- typically tclustmvmf --
%               are responsible for that check, exactly as tclust checks
%               Y once before entering the main loop rather than in
%               logmvnpdfFS).
%           mu: Mean direction. Vector. 1-by-p unit vector.
%        kappa: Concentration parameter. Scalar. kappa >= 0. kappa=0
%               corresponds to the uniform distribution on the
%               hypersphere.
%
%  Output:
%
%        logf : Log-density. Vector. n-by-1 vector containing the log of
%               the von Mises-Fisher density evaluated at each row of Y.
%
% References:
%
%   Mardia, K.V. and Jupp, P.E. (2000), "Directional Statistics", Wiley.
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   "Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions", Journal of Machine Learning Research, Vol. 6.
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('logvMFpdfFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Log-density of a circular (p=2) von Mises-Fisher distribution.
    theta = linspace(0,2*pi,50)';
    Y = [cos(theta) sin(theta)];
    mu = [1 0];
    kappa = 5;
    logf = logvMFpdfFS(Y,mu,kappa);
    plot(theta,exp(logf))
%}

%% Beginning of code

p = size(Y,2);

if kappa <= 0
    logf = repmat(vMFlogC(kappa,p), size(Y,1), 1);
    return
end

logf = vMFlogC(kappa,p) + kappa * (Y * mu');

end
%FScategory:CLUS-RobClaMULT
