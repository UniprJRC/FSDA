function logf = logvMFpdfFS(Y, mu, kappa)
%logvMFpdfFS computes the log-density of the von Mises-Fisher distribution
%
%<a href="matlab: docsearchFS('logvMFpdfFS')">Link to the help page for this function</a>
%
%   logvMFpdfFS computes, for an n-by-p matrix Y of unit vectors (rows are
%   points on the unit hypersphere S^(p-1) in R^p), the log-density of a
%   von Mises-Fisher (vMF) component with mean direction mu and
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
%   is evaluated entirely in log-space via vMFlogC, which is reliable for
%   large kappa and for large p. This is the vMF counterpart of the
%   log-space evaluation performed by logmvnpdfFS for the Gaussian density
%   used in tclust.
%
%   Several components can be evaluated in a single call: if mu is a
%   k-by-p matrix and kappa a vector of length k, the output is the
%   n-by-k matrix of the log-densities of the k components, computed with
%   one matrix product. This is how tclustvmf computes the E-step.
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
%           mu: Mean direction(s). Vector of length p, or k-by-p matrix
%               with one mean direction per row. Rows must have norm 1
%               (not enforced here; caller's responsibility).
%        kappa: Concentration parameter(s). Scalar, or vector of length k.
%               Non-negative values. kappa=0 is the uniform distribution
%               on S^(p-1). Larger kappa concentrates mass more tightly
%               around mu: the distribution is approximately
%               N(mu, (1/kappa)*I) on the tangent plane at mu for large
%               kappa.
%
%  Optional input arguments:
%
% (none beyond what is listed above; all three arguments are positional)
%
%
%  Output:
%
%        logf : Log-density values. n-by-1 vector containing
%               log f(y_i; mu, kappa) for each row y_i of Y, or n-by-k
%               matrix whose column j refers to component j when k
%               components are given. When kappa(j)<=0 all entries of
%               column j are equal to vMFlogC(0,p), the log-density of
%               the uniform distribution on S^(p-1). When kappa(j) is NaN
%               column j is NaN.
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
% See also: tclustvmf, vMFlogC, kappaFS, vMFrnd, vMFcap, sphereplot
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
    % For p=2, kappa=5, mu=[1 0], evaluate log f at 200 points on the circle
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
    %% Log-density of simulated vMF data shown on the sphere.
    % The units are divided into 4 groups according to the quartiles of
    % their log-density: the groups are concentric circular rings around
    % the mean direction.
    rng(1)
    n = 1000; mu = [1 0 0]; kappa = 20;
    Y = vMFrnd(n, mu, kappa);
    logf = logvMFpdfFS(Y, mu, kappa);
    q = quantile(logf, [0.25 0.5 0.75]);
    group = 1 + (logf<q(3)) + (logf<q(2)) + (logf<q(1));
    sphereplot(Y, 'group', group, 'labels', ...
        {'highest density','3rd quarter','2nd quarter','lowest density'});
%}
%
%{
    %% Several components in one call: discriminant analysis.
    % Units generated from two vMF distributions are assigned to the
    % component with the largest log-density. mu is a 2-by-3 matrix and
    % kappa a vector of length 2; the output has two columns.
    rng(2)
    mu = [1 0 0; 0 1 0]; kappa = [30; 10];
    Y = [vMFrnd(300, mu(1,:), kappa(1)); vMFrnd(300, mu(2,:), kappa(2))];
    logf = logvMFpdfFS(Y, mu, kappa);
    [~, group] = max(logf, [], 2);
    truegroup = [ones(300,1); 2*ones(300,1)];
    fprintf('Proportion of units correctly assigned: %.3f\n', mean(group==truegroup));
    sphereplot(Y, 'group', group, 'cap', struct('mu', mu, 'kappa', kappa));
%}
%
%{
    % Spherical (p=3) vMF: compare different concentrations.
    % The log-density is evaluated at 1000 units uniformly distributed on
    % the sphere. A large kappa gives a tightly peaked density; kappa=0
    % is the uniform density.
    rng(1)
    n = 1000;
    mu = [0 0 1];                  % north pole
    Y = vMFrnd(n, mu, 0);          % uniform on S^2
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
        mu = eye(1, p);    % arbitrary unit vector
        Y = vMFrnd(50, mu, 0);
        logf = logvMFpdfFS(Y, mu, 0);
        omega = 2 * pi^(p/2) / gamma(p/2);
        expected = -log(omega);
        fprintf('p=%2d: logf (all equal) = %.6f, expected = %.6f\n', ...
            p, logf(1), expected);
    end
%}

%% Beginning of code

[n, p] = size(Y);

% one component: mu can be a row or a column vector
if isscalar(kappa) && isvector(mu)
    mu = mu(:)';
end
kappa = kappa(:)';
k = numel(kappa);
if size(mu,1) ~= k || size(mu,2) ~= p
    error('FSDA:logvMFpdfFS:WrongInput', ...
        'mu must be a k-by-p matrix, where k is the number of elements of kappa and p the number of columns of Y')
end

% log C_p(kappa_j) in all the rows of column j (uniform density if
% kappa_j<=0)
logf = repmat(vMFlogC(kappa, p), n, 1);

% kappa_j * mu_j' * y_i for the components with positive concentration
pos = kappa > 0;
if any(pos)
    logf(:, pos) = logf(:, pos) + (Y*mu(pos,:)').*kappa(pos);
end

% components which are not defined
logf(:, isnan(kappa)) = NaN;

end
%FScategory:CLUS-RobClaMULT
