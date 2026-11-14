function [gamma1, gamma2, gamma3, kappa, beta, r, dq] = KentMstepFS(Y, w, kappa0, beta0, maxiter, tol)
%KentMstepFS estimates the parameters of a Kent distribution from a (weighted) sample
%
%<a href="matlab: docsearchFS('KentMstepFS')">Link to the help page for this function</a>
%
%   KentMstepFS estimates the parameters of a Kent(kappa,beta,Gamma)
%   distribution from a weighted sample of unit vectors on S^2. It is
%   called once per group at each M-step of the tclustkent concentration
%   loop, and it plays the role that kappaFS (together with the resultant
%   vector) plays in tclustvmf.
%
%   The estimation proceeds in two stages (Kent, 1982, Section 4;
%   Mardia et al., 2021, Section 2.3):
%
%   STAGE 1 - ORIENTATION (closed form, moment estimator of Kent 1982).
%   Let xbar and S be the weighted mean vector and the weighted matrix of
%   second moments,
%       xbar = sum_i w_i*y_i / n_w,     S = sum_i w_i*y_i*y_i' / n_w,
%   with n_w = sum_i w_i. The pole is gamma3 = xbar/||xbar||. The major
%   and minor axes gamma1 and gamma2 are the eigenvectors (largest and
%   smallest eigenvalue) of S restricted to the plane orthogonal to
%   gamma3, i.e. of the 2-by-2 matrix E'*S*E where the columns of E are
%   an orthonormal basis of that plane. In this way [gamma1 gamma2 gamma3]
%   is exactly orthogonal, gamma3'*xbar = ||xbar|| and gamma1'*S*gamma2=0.
%
%   STAGE 2 - CONCENTRATION AND OVALNESS (numerical).
%   Given the orientation, the log-likelihood divided by n_w is
%
%       kappa*r + beta*dq - log c(kappa,beta)
%
%   where r = ||xbar|| is the mean resultant length and
%   dq = gamma1'*S*gamma1 - gamma2'*S*gamma2 >= 0. This function of two
%   variables is concave and it is maximised exactly (the exact
%   normalizing constant KentlogC is used, not its normal approximation)
%   over the region kappa>=0, 0<=beta<kappa/2 with fmincon (Optimization
%   Toolbox) or, if fmincon is not available, with fminsearch.
%
%   NOTE: the Kent distribution is defined only on S^2: Y must have 3
%   columns.
%
%  Required input arguments:
%
%            Y: Input data. Matrix. n-by-3 matrix of unit vectors; each
%               row is a point on S^2. Rows are not renormalised.
%            w: Weights. Vector. n-by-1 vector of non-negative weights.
%               For crisp (CEM) assignment w is the 0/1 indicator of the
%               group; for soft (EM) assignment w contains the posterior
%               probabilities of the group. If w is empty all weights are
%               equal to 1.
%
%  Optional input arguments:
%
%       kappa0: Starting value for kappa. Scalar. Default is the vMF
%               approximation r*(3-r^2)/(1-r^2) (Banerjee et al., 2005).
%                 Example - 20
%                 Data Types - double
%        beta0: Starting value for beta. Scalar. Default is kappa0/4.
%                 Example - 5
%                 Data Types - double
%      maxiter: Maximum number of iterations of the numerical optimiser.
%               Scalar. Default 200.
%                 Example - 500
%                 Data Types - double
%          tol: Convergence tolerance of the optimiser. Scalar.
%               Default 1e-8.
%                 Example - 1e-10
%                 Data Types - double
%
%  Output:
%
%       gamma1: Major axis. Vector. 1-by-3 unit vector with the estimated
%               major axis.
%       gamma2: Minor axis. Vector. 1-by-3 unit vector with the estimated
%               minor axis.
%       gamma3: Mean direction. Vector. 1-by-3 unit vector with the
%               estimated mean direction (pole).
%        kappa: Concentration. Scalar. Estimated concentration, in
%               [0, 1e6].
%         beta: Ovalness. Scalar. Estimated ovalness, in [0, kappa/2).
%            r: Mean resultant length. Scalar. Weighted mean of gamma3'*y,
%               that is the norm of the weighted mean of the units.
%           dq: Second order statistic. Scalar. Weighted mean of
%               (gamma1'*y)^2-(gamma2'*y)^2. Together with r it forms the
%               pair of sufficient statistics of (kappa,beta) for the
%               estimated orientation, and it is needed by restrKentFS.
%
% More About:
%
% EMPTY OR DEGENERATE GROUPS
%
% When the sum of the weights is 0, or the resultant vector is null, all
% outputs are NaN: this is how an empty/degenerate group is signalled to
% tclustkent.
%
% CIRCULAR GROUPS
%
% When the two eigenvalues of the 2-by-2 matrix of Stage 1 coincide
% (dq=0) the axes gamma1 and gamma2 are not identified and the estimate
% of beta is 0: the Kent distribution reduces to the vMF distribution and
% kappa is computed from r with the vMF approximation.
%
% UPPER BOUND ON KAPPA
%
% Estimates of kappa are bounded above by 1e6. A group of (almost)
% coincident points, like some of the elemental subsets used to
% initialise tclustkent, would otherwise produce an arbitrarily large
% kappa. With kappa=1e6 the angular standard deviation of the group is
% about 0.06 degrees.
%
% References:
%
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   Clustering on the Unit Hypersphere using von Mises-Fisher Distributions,
%   J. Machine Learning Research, Vol. 6, pp. 1345-1382.
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
%   Mardia, K.V., Barber, S., Burdett, P.M., Kent, J.T. and Hamelryck, T.
%   (2021), Mixture models for spherical data with applications to protein
%   bioinformatics, arXiv:2104.13140.
%
% See also: tclustkent, KentlogC, logKentpdfFS, restrKentFS, Kentrnd, Kentcap, sphereplot
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('KentMstepFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% Recover the parameters of a Kent distribution.
    % Data are generated from a Kent distribution with pole [0 1 0],
    % elongated along [1 0 0] (an elliptical cap).
    rng(1)
    n = 2000; kappa_true = 50; beta_true = 18;
    Y = Kentrnd(n, kappa_true, beta_true, [0 1 0], [1 0 0]);
    [g1,g2,g3,kap,bet] = KentMstepFS(Y, ones(n,1));
    fprintf('kappa: true=%g, est=%.2f\n', kappa_true, kap);
    fprintf('beta:  true=%g, est=%.2f\n', beta_true, bet);
    disp('Estimated pole (true [0 1 0])'); disp(g3)
    disp('Estimated major axis (true [1 0 0], up to sign)'); disp(g1)
%}
%
%{
    %% Fitted region superimposed to the data.
    % The 95 per cent region of the fitted distribution is computed with
    % Kentcap; the units outside the region are highlighted.
    rng(2)
    n = 500;
    Y = Kentrnd(n, 80, 30, [1 0 0], [0 1 1]);
    [g1,g2,g3,kap,bet] = KentMstepFS(Y, ones(n,1));
    [P, thresh] = Kentcap(g1, g2, g3, kap, bet, 0.95);
    inside = logKentpdfFS(Y,g1,g2,g3,kap,bet) >= thresh;
    sphereplot(Y, 'group', 2-inside, 'cap', {P}, 'center', g3, ...
        'labels', {'inside','outside'});
    title(sprintf('%.1f%% of the units inside the fitted 95%% region',100*mean(inside)))
%}
%
%{
    % Weighted estimation.
    % Two elongated groups; the weights (0/1 indicators of the first
    % group) select the units used in the estimation, as in the crisp
    % M-step of tclustkent.
    rng(3)
    Y1 = Kentrnd(300, 100, 40, [1 0 0], [0 0 1]);
    Y2 = Kentrnd(300,  50, 15, [0 1 0], [1 0 0]);
    Y = [Y1; Y2];
    w = [ones(300,1); zeros(300,1)];
    [g1,g2,g3,kap,bet] = KentMstepFS(Y, w);
    fprintf('Group 1: kappa=%.1f (true 100), beta=%.1f (true 40)\n', kap, bet);
%}
%
%{
    % Circular group (beta=0, i.e. von Mises-Fisher data).
    % The estimate of beta is close to 0.
    rng(4)
    n = 2000;
    Y = Kentrnd(n, 30, 0, [0 0 1]);
    [g1,g2,g3,kap,bet] = KentMstepFS(Y, ones(n,1));
    fprintf('Circular group: kappa=%.2f (true 30), beta=%.3f (true 0)\n', kap, bet);
%}
%
%{
    % Degenerate case: empty group (all weights zero) -> NaN output.
    Y = Kentrnd(50, 20, 5);
    [g1,g2,g3,kap,bet] = KentMstepFS(Y, zeros(50,1));
    fprintf('Empty group: kappa=%.4f (should be NaN)\n', kap);
%}

%% Beginning of code

[n, p] = size(Y);
if p ~= 3
    error('FSDA:KentMstepFS:WrongDim', ...
        'KentMstepFS requires 3-column data (Kent distribution is defined on S^2 only)')
end

if nargin < 2 || isempty(w)
    w = ones(n,1);
end
w = w(:);

if nargin < 3, kappa0 = []; end
if nargin < 4, beta0 = [];  end
if nargin < 5 || isempty(maxiter), maxiter = 200; end
if nargin < 6 || isempty(tol),     tol = 1e-8;    end

% Upper bound for kappa (see More About)
kappamax = 1e6;

% Sum of the weights
nw = sum(w);

% Degenerate: empty group
if nw < 1e-10
    gamma1 = NaN(1,3); gamma2 = NaN(1,3); gamma3 = NaN(1,3);
    kappa = NaN; beta = NaN; r = NaN; dq = NaN;
    return
end

%% Stage 1: orientation

% Weighted mean vector and matrix of second moments
xbar = (w' * Y)' / nw;        % 3-by-1
S    = (Y .* w)' * Y / nw;    % 3-by-3

% mean resultant length
r = norm(xbar);

if r < 1e-10
    % null resultant: orientation not defined
    gamma1 = NaN(1,3); gamma2 = NaN(1,3); gamma3 = NaN(1,3);
    kappa = NaN; beta = NaN; r = NaN; dq = NaN;
    return
end

% pole
g3 = xbar / r;

% Orthonormal basis E=[e1 e2] of the plane orthogonal to g3
[~, imin] = min(abs(g3));
e = zeros(3,1); e(imin) = 1;
e1 = e - (e'*g3)*g3;
e1 = e1 / norm(e1);
e2 = cross(g3, e1);
E = [e1 e2];

% S restricted to the tangent plane and its spectral decomposition
T = E' * S * E;
T = (T + T')/2;
[V, D] = eig(T);
[lam, ord] = sort(diag(D), 'descend');
% major axis (largest eigenvalue) and minor axis; det([g1 g2 g3])=1
g1 = E * V(:, ord(1));
g1 = g1 / norm(g1);
g2 = cross(g3, g1);

% dq = g1'*S*g1 - g2'*S*g2 >= 0
dq = lam(1) - lam(2);

gamma1 = g1'; gamma2 = g2'; gamma3 = g3';

%% Stage 2: kappa and beta

% vMF approximation for kappa
rc = min(r, 1 - 1e-10);
kappavmf = min(rc * (3 - rc^2) / (1 - rc^2), kappamax);

% Circular group: beta=0
if dq < 1e-10
    kappa = kappavmf;
    beta = 0;
    return
end

% Starting values
if isempty(kappa0) || ~isfinite(kappa0)
    kappa0 = kappavmf;
end
if isempty(beta0) || ~isfinite(beta0)
    beta0 = kappa0 / 4;
end
kappa0 = min(max(kappa0, 0.1), kappamax);
x0 = [kappa0; max(min(beta0, kappa0/2 - 0.01), 0)];

% Minus the log-likelihood divided by nw
negpll = @(kb) -(kb(1)*r + kb(2)*dq - KentlogC(kb(1), kb(2)));

useOptimTbx = exist('fmincon', 'file') == 2;
if useOptimTbx
    % Linear inequality beta - kappa/2 <= -1e-6
    A   = [-0.5  1];
    b_c = -1e-6;
    lb  = [1e-6; 0];
    ub  = [kappamax; kappamax/2];
    opts = optimoptions('fmincon', 'Display','none', ...
        'MaxIterations', maxiter, 'FunctionTolerance', tol, ...
        'StepTolerance', tol, 'OptimalityTolerance', tol);
    try
        [xopt, ~, exitflag] = fmincon(negpll, x0, A, b_c, [], [], lb, ub, [], opts);
        if exitflag < 0 || any(~isfinite(xopt))
            xopt = x0;
        end
    catch
        xopt = x0;
    end
else
    % fminsearch on an unconstrained reparametrisation:
    % kappa = exp(t1) and beta = (kappa/2)/(1+exp(-t2)) in (0, kappa/2)
    tr = @(t) [min(exp(t(1)), kappamax); ...
        0.5*min(exp(t(1)), kappamax)/(1+exp(-t(2)))];
    ratio0 = min(max(2*x0(2)/x0(1), 1e-3), 1-1e-3);
    t0 = [log(x0(1)); log(ratio0/(1-ratio0))];
    opts = optimset('Display','none', 'MaxIter', 4*maxiter, ...
        'MaxFunEvals', 8*maxiter, 'TolFun', tol, 'TolX', tol);
    try
        topt = fminsearch(@(t) negpll(tr(t)), t0, opts);
        xopt = tr(topt);
    catch
        xopt = x0;
    end
end

kappa = min(max(xopt(1), 0), kappamax);
beta  = max(min(xopt(2), kappa/2*(1-1e-8)), 0);

end
%FScategory:CLUS-RobClaMULT
