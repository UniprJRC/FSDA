function [P, thresh] = Kentcap(gamma1, gamma2, gamma3, kappa, beta, conflev, npts)
%Kentcap computes the boundary of the confidence region of a Kent distribution on S^2
%
%<a href="matlab: docsearchFS('Kentcap')">Link to the help page for this function</a>
%
%   Kentcap is the Kent counterpart of vMFcap. The contours of constant
%   density of the von Mises-Fisher distribution are circles around the
%   mean direction, so that its confidence region is a spherical cap. The
%   contours of the Kent(kappa,beta,Gamma) distribution are instead ovals
%   around the pole gamma3, elongated along the major axis gamma1: the
%   confidence region is an "elliptical cap". Kentcap returns npts points
%   on the boundary of the region
%
%     {y in S^2 : f(y; kappa,beta,Gamma) >= t}
%
%   where the threshold t is such that the probability of the region
%   under Kent(kappa,beta,Gamma) is conflev. The points are ordered along
%   the contour, which is a closed curve.
%
%  Required input arguments:
%
%       gamma1: Major axis. Vector of length 3. Its component along
%               gamma3 is removed and the remainder is rescaled to unit
%               norm (as in Kentrnd).
%       gamma2: Minor axis. Vector of length 3 or empty value. If it is
%               empty, gamma2=cross(gamma3,gamma1) is used. It is
%               supplied only to have the same list of arguments as
%               logKentpdfFS: the region does not depend on its sign.
%       gamma3: Mean direction (pole). Vector of length 3.
%        kappa: Concentration parameter. Scalar. kappa>=0.
%         beta: Ovalness parameter. Scalar. 0<=beta<kappa/2. With beta=0
%               the region is the circular cap of vMFcap.
%
%  Optional input arguments:
%
%      conflev: Confidence level. Scalar in (0,1). Default 0.95.
%                 Example - 0.90
%                 Data Types - double
%         npts: Number of points on the boundary. Scalar. Default 200.
%               The first and the last point coincide (closed curve).
%                 Example - 300
%                 Data Types - double
%
%  Output:
%
%            P: Boundary points. npts-by-3 matrix of unit vectors, ordered
%               along the contour. Plot it with
%               plot3(P(:,1),P(:,2),P(:,3)) or pass it to sphereplot.
%       thresh: Value of the log-density on the boundary. Scalar. A point
%               y is inside the region if
%               logKentpdfFS(y,gamma1,gamma2,gamma3,kappa,beta)>=thresh.
%
% More About:
%
% In the canonical frame, with colatitude theta (angle from the pole) and
% longitude phi (measured from the major axis), the log-density is, up to
% the normalizing constant,
%
%   g(theta,phi) = kappa*cos(theta) + beta*sin(theta)^2*cos(2*phi).
%
% The probability of the region {g>=t} is computed by numerical
% integration on a fine grid of (cos(theta),phi), and t is the value for
% which this probability is conflev. For each longitude the equation
% g=t is quadratic in c=cos(theta),
%
%   A*c^2 - kappa*c + (t-A) = 0,     A = beta*cos(2*phi),
%
% and, because g is increasing in c when kappa>2*beta, it has a single
% root in [-1,1]. For concentrated distributions the region is
% approximately an ellipse with angular semi-axes proportional to
% 1/sqrt(kappa-2*beta) (along gamma1) and 1/sqrt(kappa+2*beta) (along
% gamma2).
%
% References:
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
% See also: vMFcap, Kentrnd, sphereplot, logKentpdfFS, tclustkent
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('Kentcap')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% Elliptical caps at different confidence levels.
    rng(1)
    g3 = [1 0 0]; g1 = [0 0 1]; kappa = 60; beta = 25;
    Y = Kentrnd(500, kappa, beta, g3, g1);
    P50 = Kentcap(g1, [], g3, kappa, beta, 0.50);
    P95 = Kentcap(g1, [], g3, kappa, beta, 0.95);
    P99 = Kentcap(g1, [], g3, kappa, beta, 0.99);
    % the three regions refer to the same group (same colour)
    sphereplot(Y, 'cap', {P50, P95, P99}, 'center', [g3; g3; g3], ...
        'capgroup', [1 1 1], 'capfill', false);
%}
%
%{
    % Coverage of the region.
    % The proportion of simulated units which fall inside the region is
    % close to the confidence level.
    rng(2)
    g3 = [0 1 0]; g1 = [1 0 0]; kappa = 40; beta = 15;
    [Y, Gamma] = Kentrnd(100000, kappa, beta, g3, g1);
    for conflev = [0.5 0.9 0.95 0.99]
        [~, thresh] = Kentcap(g1, [], g3, kappa, beta, conflev);
        logf = logKentpdfFS(Y, Gamma(:,1), Gamma(:,2), Gamma(:,3), kappa, beta);
        fprintf('conflev=%.2f: proportion inside %.4f\n', conflev, mean(logf>=thresh));
    end
%}
%
%{
    % With beta=0 the region is the circular cap of vMFcap.
    mu = [0 0 1]; kappa = 20;
    P1 = Kentcap([1 0 0], [], mu, kappa, 0);
    P2 = vMFcap(mu, kappa);
    fprintf('Angular radius: Kentcap %.4f, vMFcap %.4f (degrees)\n', ...
        acosd(P1(1,:)*mu'), acosd(P2(1,:)*mu'));
%}

%% Beginning of code

if nargin < 5
    error('FSDA:Kentcap:missingInputs', 'gamma1, gamma2, gamma3, kappa and beta must be supplied')
end
if nargin < 6 || isempty(conflev)
    conflev = 0.95;
end
if nargin < 7 || isempty(npts)
    npts = 200;
end
if conflev <= 0 || conflev >= 1
    error('FSDA:Kentcap:WrongConflev', 'conflev must be a scalar in the interval (0,1)')
end

% Orthonormal frame (g1 major axis, g2 minor axis, g3 pole)
g3 = gamma3(:)';
g3 = g3/norm(g3);
g1 = gamma1(:)';
g1 = g1 - (g1*g3')*g3;
if norm(g1) < 1e-10
    error('FSDA:Kentcap:WrongAxis', 'gamma1 must not be parallel to gamma3')
end
g1 = g1/norm(g1);
if isempty(gamma2)
    g2 = cross(g3, g1);
else
    g2 = gamma2(:)';
    g2 = g2 - (g2*g3')*g3 - (g2*g1')*g1;
    if norm(g2) < 1e-10
        g2 = cross(g3, g1);
    else
        g2 = g2/norm(g2);
    end
end

phi = linspace(0, 2*pi, npts)';

if kappa <= 0
    % uniform distribution: cap of area conflev around the pole
    cb = repmat(1 - 2*conflev, npts, 1);
    thresh = -log(4*pi);
else
    % 1) Threshold: numerical integration on a grid of d=1-cos(theta) and
    % phi (the area element is dd*dphi). The density is negligible (less
    % than exp(-40) times its maximum) for d>40/(kappa-2*beta).
    nd = 2000; nph = 360;
    dmax = min(2, 40/max(kappa-2*beta, 1e-8*kappa));
    d = ((1:nd)' - 0.5)*dmax/nd;
    ph = ((1:nph) - 0.5)*2*pi/nph;
    c = 1 - d;
    g = kappa*c + beta*(1 - c.^2).*cos(2*ph);     % nd-by-nph
    gs = sort(g(:), 'descend');
    cumw = cumsum(exp(gs - kappa));
    cumw = cumw/cumw(end);
    t = gs(find(cumw >= conflev, 1));

    % 2) Boundary: root in [-1,1] of A*c^2 - kappa*c + (t-A) = 0
    A = beta*cos(2*phi);
    cb = repmat(t/kappa, npts, 1);
    nz = abs(A) > 1e-12*kappa;
    cb(nz) = (kappa - sqrt(max(kappa^2 - 4*A(nz).*(t - A(nz)), 0)))./(2*A(nz));
    cb = min(max(cb, -1), 1);
    thresh = t - KentlogC(kappa, beta);
end

sb = sqrt(1 - cb.^2);
P = cb*g3 + (sb.*cos(phi))*g1 + (sb.*sin(phi))*g2;

end
%FScategory:CLUS-RobClaMULT
