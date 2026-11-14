function logf = logKentpdfFS(Y, gamma1, gamma2, gamma3, kappa, beta)
%logKentpdfFS computes the log-density of the Kent distribution on S^2
%
%<a href="matlab: docsearchFS('logKentpdfFS')">Link to the help page for this function</a>
%
%   logKentpdfFS evaluates, for an n-by-3 matrix Y of unit vectors on S^2,
%   the log-density of a single Kent (Fisher-Bingham FB5) component:
%
%       log f(x; kappa, beta, Gamma) = -log c(kappa,beta)
%           + kappa*(gamma3'*x) + beta*[(gamma1'*x)^2 - (gamma2'*x)^2]
%
%   where:
%     gamma3 is the 3-by-1 mean direction (pole) unit vector,
%     gamma1 is the 3-by-1 major-axis unit vector (gamma1 perp gamma3),
%     gamma2 = cross(gamma3,gamma1) is the minor-axis unit vector,
%     kappa  >= 0 is the concentration parameter,
%     beta   >= 0 is the ovalness parameter (must satisfy beta < kappa/2),
%     c(kappa,beta) is the normalizing constant computed by KentlogC.
%
%   The notation follows Kent (1982, eq. 2) and Mardia et al. (2021, eq. 2),
%   but with the axis labelling convention of Dong & Melnykov (2024):
%     gamma1 = major axis  (direction of maximum secondary concentration)
%     gamma2 = minor axis  (direction of minimum secondary concentration)
%     gamma3 = mean direction / pole
%
%   This is the Kent-distribution counterpart of logvMFpdfFS: it plays the
%   same role in tclustkent that logvMFpdfFS plays in tclustvmf — it
%   provides the per-unit log-discriminant values that enter the E-step.
%
%   The special case beta=0 reduces to the vMF density with parameter kappa:
%   the quadratic term vanishes and the formula coincides with logvMFpdfFS
%   (up to a constant difference from the different normalizing constants,
%   which is absorbed into logc(kappa,0); see KentlogC).
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-3 matrix of unit vectors; each row
%               is a point on S^2. Rows are NOT renormalized inside this
%               function: the caller (tclustkent) is responsible for the
%               unit-norm check.
%       gamma1: Major-axis direction. Vector. 3-by-1 (or 1-by-3) unit
%               vector; must be orthogonal to gamma3. Caller's responsibility.
%       gamma2: Minor-axis direction. Vector. 3-by-1 (or 1-by-3) unit
%               vector; must be orthogonal to both gamma1 and gamma3.
%               Typically gamma2 = cross(gamma3, gamma1).
%       gamma3: Mean direction (pole). Vector. 3-by-1 (or 1-by-3) unit
%               vector; the mode of the density.
%        kappa: Concentration parameter. Scalar. Non-negative.
%               kappa=0 is the uniform distribution on S^2.
%         beta: Ovalness parameter. Scalar. Non-negative; must satisfy
%               0 <= beta < kappa/2 for a proper density. beta=0 gives
%               the isotropic (vMF) special case.
%
%  Optional input arguments:
%
% (none beyond what is listed above; all arguments are positional)
%
%  Output:
%
%        logf : Log-density values. Vector. Column vector of length n with
%               the values log f(y_i; kappa, beta, Gamma) for each row
%               y_i of Y. When kappa <= 0, all entries equal
%               log(1/(4*pi)) (uniform).
%
% More About:
%
% ORIENTATION CONVENTION
%
% The 3-by-3 orientation matrix is Gamma = [gamma1 gamma2 gamma3] where the
% columns are the major axis, minor axis, and mean direction (pole). This
% matches Dong & Melnykov (2024), where gamma1 is the mean and gamma1,
% gamma2, gamma3 are ordered as (mean, major, minor). We follow the
% Mardia et al. (2021) / Kent (1982) convention where gamma3 is the mean
% direction (pole) and gamma1, gamma2 are the ovalness axes, since this
% makes the quadratic form diagonal in the rotated frame.
%
% PARAMETER COUNT
%
% A Kent distribution on S^2 has 5 free parameters:
%   kappa (1), beta (1), and the orientation matrix Gamma (3 rotation
%   angles for SO(3), minus 1 internal symmetry of the ovalness axes,
%   = 3; total = 5). This agrees with the name "FB5" (Fisher-Bingham
%   5-parameter distribution). The free parameters of Gamma are
%   typically parameterised as (phi, eta, psi) in Kent (1982).
%
% References:
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
%   Mardia, K.V., Barber, S., Burdett, P.M., Kent, J.T. and Hamelryck, T.
%   (2021), Mixture models for spherical data with applications to protein
%   bioinformatics, arXiv:2104.13140.
%
%   Dong, A. and Melnykov, V. (2024), Contaminated Kent mixture model for
%   clustering non-spherical directional data with heavy tails or scatter,
%   Statistics and Probability Letters, Vol. 208, 110058.
%
% See also: tclustkent, KentlogC, KentMstepFS, restrKentFS, Kentrnd, Kentcap, sphereplot
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('logKentpdfFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Log-density at the mode: at x=gamma3 the quadratic term is 0 and
    % the vMF term is kappa. So logf(gamma3) = -logc + kappa.
    kappa = 20; beta = 5;
    gamma3 = [0 0 1]; gamma1 = [1 0 0]; gamma2 = [0 1 0];
    logf_mode = logKentpdfFS(gamma3, gamma1, gamma2, gamma3, kappa, beta);
    logf_expected = -KentlogC(kappa,beta) + kappa;
    fprintf('logf at mode = %.6f, expected = %.6f\n', logf_mode, logf_expected);
%}
%
%{
    % For beta=0 the Kent density reduces to the vMF density.
    % Check that logKentpdfFS and logvMFpdfFS agree (after accounting for
    % the different normalizing constants).
    rng(1)
    n = 100; kappa = 15; beta = 0;
    % points uniformly distributed on the sphere (vMF with kappa=0)
    Y = vMFrnd(n, [0 0 1], 0);
    gamma3 = [0 0 1]; gamma1 = [1 0 0]; gamma2 = [0 1 0];
    logf_kent = logKentpdfFS(Y, gamma1, gamma2, gamma3, kappa, beta);
    logf_vmf  = logvMFpdfFS(Y, gamma3, kappa);
    % should agree to machine precision
    fprintf('max|logf_kent - logf_vmf| = %.2e\n', max(abs(logf_kent - logf_vmf)));
%}
%
%{
    % Density is highest near the pole (gamma3) and varies with beta.
    % With larger beta the density becomes more elongated along gamma1.
    rng(2)
    n = 500; kappa = 30;
    % points uniformly distributed on the sphere (vMF with kappa=0)
    Y = vMFrnd(n, [0 0 1], 0);
    gamma3 = [0 0 1]; gamma1 = [1 0 0]; gamma2 = [0 1 0];
    for beta = [0 5 10 14]
        logf = logKentpdfFS(Y, gamma1, gamma2, gamma3, kappa, beta);
        fprintf('beta=%2d: mean logf=%.3f, range=[%.3f %.3f]\n', ...
            beta, mean(logf), min(logf), max(logf));
    end
%}
%
%{
    %% Log-density of simulated Kent data shown on the sphere.
    % The units are divided into 4 groups according to the quartiles of
    % their log-density: the groups are concentric elliptical rings.
    rng(3)
    n = 1000; kappa = 40; beta = 15;
    [Y, Gamma] = Kentrnd(n, kappa, beta, [1 0 0], [0 0 1]);
    logf = logKentpdfFS(Y, Gamma(:,1), Gamma(:,2), Gamma(:,3), kappa, beta);
    q = quantile(logf, [0.25 0.5 0.75]);
    group = 1 + (logf<q(3)) + (logf<q(2)) + (logf<q(1));
    sphereplot(Y, 'group', group, 'labels', ...
        {'highest density','3rd quarter','2nd quarter','lowest density'});
%}

%% Beginning of code

n = size(Y, 1);

% Ensure column vectors for inner products
g1 = gamma1(:);
g2 = gamma2(:);
g3 = gamma3(:);

if kappa <= 0
    % Uniform distribution on S^2
    logf = repmat(-log(4*pi), n, 1);
    return
end

logc = KentlogC(kappa, beta);

% Log-density: -logc + kappa*(gamma3'*x) + beta*[(gamma1'*x)^2 - (gamma2'*x)^2]
dot3 = Y * g3;          % n-by-1: gamma3' * x_i
dot1 = Y * g1;          % n-by-1: gamma1' * x_i
dot2 = Y * g2;          % n-by-1: gamma2' * x_i

logf = -logc + kappa * dot3 + beta * (dot1.^2 - dot2.^2);

end
%FScategory:CLUS-RobClaMULT
