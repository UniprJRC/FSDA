function [Y, Gamma] = Kentrnd(n, kappa, beta, gamma3, gamma1)
%Kentrnd generates random unit vectors from the Kent distribution on S^2
%
%<a href="matlab: docsearchFS('Kentrnd')">Link to the help page for this function</a>
%
%   Kentrnd returns an n-by-3 matrix whose rows are independent draws from
%   the Kent (Fisher-Bingham FB5) distribution Kent(kappa,beta,Gamma) on
%   the unit sphere S^2, with density
%
%     f(x) = (1/c(kappa,beta)) *
%            exp{ kappa*(gamma3'*x) + beta*[(gamma1'*x)^2 - (gamma2'*x)^2] }
%
%   where gamma3 is the mean direction (pole), gamma1 the major axis and
%   gamma2 = cross(gamma3,gamma1) the minor axis. The sample is exact (no
%   approximation is involved): it is obtained by acceptance-rejection
%   from a von Mises-Fisher envelope (see "More About").
%
%   For beta=0 the draws are von Mises-Fisher with concentration kappa
%   (circular contours around gamma3). For beta>0 the contours of constant
%   density are ovals elongated along gamma1: the data form an "elliptical
%   cap", whose axis ratio is approximately
%   sqrt((kappa+2*beta)/(kappa-2*beta)) for concentrated data.
%
%  Required input arguments:
%
%            n: Number of random vectors. Scalar. Positive integer.
%        kappa: Concentration parameter. Scalar. kappa>=0.
%         beta: Ovalness parameter. Scalar. 0<=beta<kappa/2 (the condition
%               for a unimodal Kent density). beta=0 gives the von
%               Mises-Fisher distribution.
%
%  Optional input arguments:
%
%       gamma3: Mean direction (pole). Vector of length 3. It is rescaled
%               to unit norm. Default is [0 0 1].
%                 Example - [1 0 0]
%                 Data Types - double
%       gamma1: Major axis. Vector of length 3. Direction along which the
%               cluster is elongated. It does not need to be exactly
%               orthogonal to gamma3: its component along gamma3 is
%               removed and the remainder is rescaled to unit norm. It
%               must not be parallel to gamma3. Default is an arbitrary
%               direction orthogonal to gamma3 (irrelevant when beta=0).
%                 Example - [0 0 1]
%                 Data Types - double
%
%  Output:
%
%            Y: Random sample. Matrix. n-by-3 matrix of unit vectors, one
%               draw from Kent(kappa,beta,Gamma) per row.
%        Gamma: Orientation matrix actually used. Matrix. 3-by-3 orthogonal
%               matrix [gamma1 gamma2 gamma3] (columns: major axis, minor
%               axis, pole), after the orthonormalisation described above.
%
% More About:
%
% ACCEPTANCE-REJECTION FROM A VON MISES-FISHER ENVELOPE
%
% In the canonical frame (pole = third axis) write x=(x1,x2,x3). The
% proposal is the vMF distribution with the same pole and concentration
% a = kappa-2*beta > 0, which on S^2 is sampled exactly by inversion:
%
%     x3 = 1 + log( u + (1-u)*exp(-2*a) )/a,    u ~ U(0,1),
%
% with a uniform longitude. A proposed point is accepted with probability
%
%     exp{ -2*beta*(1-x3) + beta*(x1^2-x2^2) }.
%
% This is a valid probability because x1^2-x2^2 <= 1-x3^2, so that the
% exponent is at most -beta*(1-x3)^2 <= 0; and it is the ratio between the
% Kent kernel and exp(2*beta) times the vMF kernel, so that accepted
% points are exactly Kent distributed. For concentrated data the
% acceptance rate is approximately sqrt((kappa-2*beta)/(kappa+2*beta)): it
% is 1 for beta=0 and decreases as the contours become more elongated.
%
% References:
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
%   Kent, J.T., Ganeiber, A.M. and Mardia, K.V. (2018), A new unified
%   approach for the simulation of a wide class of directional
%   distributions, Journal of Computational and Graphical Statistics,
%   Vol. 27, pp. 291-301.
%
% See also: vMFrnd, sphereplot, Kentcap, tclustkent, logKentpdfFS, KentMstepFS, KentlogC
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('Kentrnd')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% An elliptical cap on the sphere.
    % 500 draws from a Kent distribution with pole [1 0 0], elongated
    % along the direction [0 0 1].
    rng(1)
    g3 = [1 0 0]; g1 = [0 0 1]; kappa = 100; beta = 40;
    Y = Kentrnd(500, kappa, beta, g3, g1);
    % data and 95 per cent confidence region
    sphereplot(Y, 'cap', struct('gamma3',g3,'gamma1',g1,'kappa',kappa,'beta',beta));
    title('Kent(\kappa=100, \beta=40): elliptical cap')
%}
%
%{
    %% Effect of beta: from a circular to an elongated cap.
    % kappa is fixed; the cap is circular for beta=0 (von Mises-Fisher)
    % and becomes more and more elongated as beta approaches kappa/2.
    rng(2)
    kappa = 60;
    betas = [0 15 25];
    figure
    for i = 1:3
        Y = Kentrnd(400, kappa, betas(i), [1 0 0], [0 0 1]);
        sphereplot(Y, 'markersize', 8, 'h', subplot(1,3,i));
        title(sprintf('\\kappa=%g, \\beta=%g',kappa,betas(i)))
    end
%}
%
%{
    % Check of the sample moments in the canonical frame.
    % For concentrated data the coordinates along the major and minor
    % axes are approximately normal with variances 1/(kappa-2*beta) and
    % 1/(kappa+2*beta) (Kent, 1982).
    rng(3)
    kappa = 200; beta = 60;
    [Y, Gamma] = Kentrnd(100000, kappa, beta, [1 1 0], [0 0 1]);
    Z = Y*Gamma;            % coordinates in the frame (major, minor, pole)
    fprintf('var major axis: %.5f (approx. %.5f)\n', var(Z(:,1)), 1/(kappa-2*beta));
    fprintf('var minor axis: %.5f (approx. %.5f)\n', var(Z(:,2)), 1/(kappa+2*beta));
%}
%
%{
    % The parameters used for the simulation are recovered by KentMstepFS.
    rng(4)
    Y = Kentrnd(2000, 50, 18, [0 1 0], [1 0 0]);
    [g1,g2,g3,kap,bet] = KentMstepFS(Y, ones(2000,1));
    fprintf('kappa: true 50, estimated %.2f\n', kap);
    fprintf('beta:  true 18, estimated %.2f\n', bet);
    disp('Estimated pole and major axis (up to sign)'); disp([g3; g1])
%}

%% Beginning of code

if nargin < 3
    error('FSDA:Kentrnd:missingInputs', 'n, kappa and beta must be supplied')
end
if kappa < 0 || beta < 0
    error('FSDA:Kentrnd:WrongPar', 'kappa and beta must be non-negative')
end
if beta > 0 && 2*beta >= kappa
    error('FSDA:Kentrnd:WrongPar', 'beta must be smaller than kappa/2')
end

if nargin < 4 || isempty(gamma3)
    gamma3 = [0 0 1];
end
g3 = gamma3(:)';
g3 = g3/norm(g3);

if nargin < 5 || isempty(gamma1)
    [~, imin] = min(abs(g3));
    gamma1 = zeros(1,3);
    gamma1(imin) = 1;
end
g1 = gamma1(:)';
g1 = g1 - (g1*g3')*g3;
if norm(g1) < 1e-10
    error('FSDA:Kentrnd:WrongAxis', 'gamma1 must not be parallel to gamma3')
end
g1 = g1/norm(g1);
g2 = cross(g3, g1);
Gamma = [g1' g2' g3'];

% Concentration of the vMF envelope
a = kappa - 2*beta;

X = zeros(n, 3);
nacc = 0;
while nacc < n
    % Number of proposals: missing draws divided by the (approximate)
    % acceptance rate, with some slack
    m = ceil(1.2*(n-nacc)*sqrt((kappa+2*beta)/max(a,eps))) + 20;
    m = min(m, 1e6);

    % vMF(a) proposal in the canonical frame (exact inversion on S^2)
    u = rand(m,1);
    if a > 1e-8
        w = 1 + log(u + (1-u).*exp(-2*a))/a;
    else
        w = 2*u - 1;
    end
    w = min(max(w,-1),1);
    phi = 2*pi*rand(m,1);
    s = sqrt(1 - w.^2);
    x1 = s.*cos(phi);
    x2 = s.*sin(phi);

    % Accept with probability exp(-2*beta*(1-x3) + beta*(x1^2-x2^2))
    acc = log(rand(m,1)) <= -2*beta*(1-w) + beta*(x1.^2 - x2.^2);
    idx = find(acc, n-nacc);
    na = numel(idx);
    X(nacc+1:nacc+na, :) = [x1(idx) x2(idx) w(idx)];
    nacc = nacc + na;
end

% Rotate from the canonical frame to (gamma1, gamma2, gamma3)
Y = X(:,1)*g1 + X(:,2)*g2 + X(:,3)*g3;

end
%FScategory:CLUS-RobClaMULT
