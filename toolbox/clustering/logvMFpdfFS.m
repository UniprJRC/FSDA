function logf = logvMFpdfFS(Y, mu, kappa)
%logvMFpdfFS computes the log-density of the von Mises-Fisher distribution
%
%<a href="matlab: docsearchFS('logvMFpdfFS')">Link to the help page for this function</a>
%
%   logvMFpdfFS computes, for an n-by-p matrix Y of unit vectors (rows are
%   points on the unit hypersphere S^(p-1) in R^p), the log-density of a
%   single von Mises-Fisher (vMF) component with mean direction mu and
%   concentration kappa:
%
%       log f(y; mu, kappa) = log C_p(kappa) + kappa * mu' * y
%
%   where the normalizing constant is
%
%       C_p(kappa) = kappa^(p/2-1) / ( (2*pi)^(p/2) * I_(p/2-1)(kappa) )
%
%   and I_nu denotes the modified Bessel function of the first kind of
%   order nu. Special cases: p=2 gives the von Mises distribution on the
%   circle (I_{0}(kappa) in the denominator); p=3 gives the spherical vMF
%   with C_3(kappa) = kappa / (4*pi*sinh(kappa)). The normalizing constant
%   is evaluated entirely in log-space via vMFlogC, which uses the
%   exponentially-scaled Bessel function besseli(...,1) to avoid overflow
%   for large kappa or large p; the scaling factor cancels in the log. This
%   is the vMF counterpart of the log-space evaluation performed by
%   logmvnpdfFS for the Gaussian density used in tclust.
%
%   The special case kappa=0 (or kappa<0, treated as 0) returns the
%   log-density of the uniform distribution on S^(p-1), which equals
%   log(1/omega_{p-1}) = -log(omega_{p-1}) for all y, where
%   omega_{p-1} = 2*pi^(p/2)/Gamma(p/2) is the surface area of S^(p-1).
%   This is consistent with the vMF model since A_p(0)=0 (mean resultant
%   length zero means maximum dispersion) and C_p(0) = 1/omega_{p-1}.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-p matrix of unit vectors; each
%               row is a point on S^(p-1). Rows are NOT renormalized
%               inside this function: the caller (typically tclustvmf) is
%               responsible for the unit-norm check, exactly as tclust
%               checks Y once before the main loop rather than inside
%               logmvnpdfFS.
%           mu: Mean direction. Vector. 1-by-p unit vector specifying the
%               centre of the vMF distribution on S^(p-1). Must have
%               norm 1 (not enforced here; caller's responsibility).
%        kappa: Concentration parameter. Scalar. Non-negative value.
%               kappa=0 is the uniform distribution on S^(p-1). Larger
%               kappa concentrates mass more tightly around mu: the
%               distribution is approximately N(mu, (1/kappa)*I) on the
%               tangent plane at mu for large kappa.
%
%  Optional input arguments:
%
% (none beyond what is listed above; all three arguments are positional)
%
%
%  Output: 
%
%        logf : Log-density values. Vector. n-by-1 vector containing
%               log f(y_i; mu, kappa) for each row y_i of Y. When
%               kappa<=0 all entries are equal to vMFlogC(0,p), the
%               log-density of the uniform distribution on S^(p-1).
%
% References:
%
%   Mardia, K.V. and Jupp, P.E. (2000), Directional Statistics, Wiley,
%   Chichester. (vMF density and normalizing constant)
%
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions, Journal of Machine Learning Research, Vol. 6,
%   pp. 1345-1382. (log-space evaluation of C_p(kappa))
%
% See also: tclustvmf, vMFlogC, kappaFS
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('logvMFpdfFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Log-density of a circular (p=2) vMF distribution (von Mises).
    % For p=2, kappa=5, mu=[1 0], evaluate log f at 50 points on the circle
    % and plot the density. The mode is at theta=0 (direction of mu).
    theta = linspace(0, 2*pi, 200)';
    Y = [cos(theta), sin(theta)];
    mu = [1 0];
    kappa = 5;
    logf = logvMFpdfFS(Y, mu, kappa);
    figure;
    plot(theta, exp(logf), 'b-', 'LineWidth', 1.5);
    xlabel('theta'); ylabel('density');
    title(sprintf('von Mises density, kappa=%.0f', kappa));
%}
%
%{
    % Spherical (p=3) vMF: compare two concentrations.
    % A large kappa gives a tightly peaked density; kappa=0 is uniform.
    % Evaluate log f at 1000 random unit vectors on S^2.
    rng(1)
    n = 1000;
    Y = randn(n, 3);
    Y = Y ./ vecnorm(Y, 2, 2);    % project to S^2
    mu = [0 0 1];                  % north pole
    for kappa = [0, 1, 10, 50]
        logf = logvMFpdfFS(Y, mu, kappa);
        fprintf('kappa=%3d: mean logf = %8.3f, range = [%.3f, %.3f]\n', ...
            kappa, mean(logf), min(logf), max(logf));
    end
%}
%
%{
    % Verify the kappa=0 case returns the uniform log-density.
    % For any p, f(y; mu, 0) = 1/omega_{p-1} for all y,
    % where omega_{p-1} = 2*pi^(p/2)/Gamma(p/2) is the surface area of S^(p-1).
    rng(1)
    for p = [2, 3, 5, 10]
        Y = randn(50, p);
        Y = Y ./ vecnorm(Y, 2, 2);
        mu = eye(1, p);    % arbitrary unit vector
        logf = logvMFpdfFS(Y, mu, 0);
        omega = 2 * pi^(p/2) / gamma(p/2);
        expected = -log(omega);
        fprintf('p=%2d: logf (all equal) = %.6f, expected = %.6f\n', ...
            p, logf(1), expected);
    end
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
