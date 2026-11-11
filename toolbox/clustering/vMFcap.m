function P = vMFcap(mu,kappa,conflev,npts)
%vMFcap computes the boundary of a von Mises-Fisher confidence cap in any dimension
%
%<a href="matlab: docsearchFS('vMFcap')">Link to the help page for this function</a>
%
%   vMFcap is the directional-data counterpart of a Gaussian confidence
%   ellipse: because the vMF(mu,kappa) density depends on y only through
%   the inner product mu'*y (equivalently, the angle between y and mu),
%   its iso-density contours are exactly the spherical caps
%   {y in S^(p-1) : angle(y,mu) <= theta}, for any embedding dimension p.
%   vMFcap returns npts points on the boundary of the cap whose probability
%   under vMF(mu,kappa) equals conflev, i.e. it returns the boundary circle
%   (or its Monte Carlo approximation for p~=3) of
%
%     {y : mu'*y >= t_alpha}   where P(mu'*y >= t_alpha) = conflev.
%
%   For p=3 the boundary is a genuine circle (a 1-sphere embedded in R^3)
%   and vMFcap returns it exactly, as npts evenly spaced points. For p~=3
%   the boundary is a (p-2)-sphere in the (p-1)-dimensional subspace
%   orthogonal to mu; vMFcap returns npts points drawn uniformly on that
%   submanifold (in the native R^p) as a Monte Carlo representation. No
%   dimensionality reduction is performed here: how to view a
%   high-dimensional cap in 2 or 3 dimensions is left to the caller (see
%   tclustvmf's 'plotdims' option).
%
%  Required input arguments:
%
%          mu: Mean direction. Vector. 1-by-p unit vector; the centre of
%             the vMF distribution whose confidence cap is requested.
%       kappa: Concentration parameter. Scalar. kappa>=0. kappa=0 gives
%             the uniform distribution on S^(p-1); the resulting cap is
%             then a spherical cap of opening angle acos(1-2*conflev).
%
%  Optional input arguments:
%
%     conflev: Confidence level of the cap. Scalar in (0,1). Default 0.95.
%              The cap boundary is the set of directions whose angular
%              distance from mu equals acos(t_alpha), where t_alpha is the
%              (1-conflev) quantile of the marginal distribution of mu'*y
%              under vMF(mu,kappa). See "More About" for the formula.
%                 Example - 0.90
%                 Data Types - double
%       npts : Number of points used to represent the boundary. Scalar.
%              Default 100. For p==3 these are evenly spaced around the
%              exact circle. For p~=3 they are a random uniform sample of
%              the boundary (p-2)-sphere in the subspace orthogonal to mu.
%                 Example - 200
%                 Data Types - double
%
%  Output:
%
%           P: Boundary points. Matrix. npts-by-p matrix of unit vectors
%              on the boundary of the confidence cap. For p==3 all rows
%              are exactly on the circle in S^2; plot with
%              plot3(P(:,1),P(:,2),P(:,3)). For p~=3 P is a Monte Carlo
%              sample of the boundary (p-2)-sphere; project P (and the
%              data) onto a 3D view before plotting.
%
% More About:
%
% ANGULAR RADIUS OF THE CAP
%
% The marginal distribution of t = mu'*y under vMF(mu,kappa) has density
% proportional to (Mardia and Jupp, 2000, Theorem 9.3.5):
%
%   g(t) = exp(kappa*t) * (1 - t^2)^((p-3)/2),   t in [-1,1].
%
% The cap threshold t_alpha satisfies P(t >= t_alpha) = conflev, i.e.
% t_alpha is the (1-conflev) quantile of g. The angular radius is then
% theta_alpha = acos(t_alpha).
%
% For p=3 the exponent (p-3)/2 = 0 and g reduces to exp(kappa*t), giving
% the closed form (rearranged for numerical stability at large kappa):
%
%   t_alpha = 1 + log( (1-conflev) + conflev*exp(-2*kappa) ) / kappa.
%
% As kappa -> inf, exp(-2*kappa) -> 0 so the argument of log approaches
% 1-conflev > 0, keeping t_alpha well defined. For kappa=0 the formula
% reduces to 1 + log(1)/0, which is 0/0, so the kappa=0 branch is
% handled separately (see below).
%
% For p~=3 there is no closed form, so t_alpha is found by numerical
% integration of g (with the stabilised form exp(kappa*(t-1))*(1-t^2)^((p-3)/2)
% that keeps the integrand in [0,1]) followed by a 1-D root-finding step.
%
% For kappa=0 (uniform on S^(p-1), any p), t = mu'*y is uniform on [-1,1]
% regardless of p, so t_alpha = 1 - 2*conflev exactly.
%
% BOUNDARY PARAMETRISATION
%
% Once t_alpha = cos(theta_alpha) is known, every point on the boundary is
%
%   P_i = cos(theta_alpha)*mu + sin(theta_alpha)*v_i
%
% where v_i is a unit vector in the (p-1)-dimensional subspace orthogonal
% to mu. This satisfies ||P_i|| = 1 (unit sphere) and mu'*P_i = cos(theta_alpha)
% = t_alpha for all i (exactly on the cap boundary). For p=3 the subspace
% is 2-dimensional and v_i = cos(phi_i)*e1 + sin(phi_i)*e2 traces the
% exact circle as phi_i sweeps [0, 2*pi]. For general p, v_i is drawn
% uniformly from the (p-2)-sphere by projecting Gaussian noise off mu.
%
% References:
%
%   Mardia, K.V. and Jupp, P.E. (2000), Directional Statistics, Wiley,
%   Chichester. (Theorem 9.3.5: marginal distribution of mu'*y under vMF)
%
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions, Journal of Machine Learning Research, Vol. 6,
%   pp. 1345-1382.
%
% See also: tclustvmf, logvMFpdfFS, kappaFS
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('vMFcap')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % 95% confidence caps for two vMF groups on S^2 (p=3).
    % A tight group (kappa=40) around [1 0 0] and a diffuse group
    % (kappa=8) around [0 1 0]. The caps are drawn as exact circles.
    mu1=[1 0 0]; kappa1=40;
    mu2=[0 1 0]; kappa2=8;
    P1 = vMFcap(mu1,kappa1);
    P2 = vMFcap(mu2,kappa2);
    [Yx,Yy,Yz] = sphere(40);
    figure;
    surf(Yx,Yy,Yz,'FaceColor',[0.9 0.9 1],'EdgeColor','none','FaceAlpha',0.25);
    hold on
    plot3(P1(:,1),P1(:,2),P1(:,3),'b-','LineWidth',2)
    plot3(P2(:,1),P2(:,2),P2(:,3),'r-','LineWidth',2)
    axis equal; view(3)
    legend({'','kappa=40','kappa=8'})
%}
%
%{
    %% Effect of confidence level on the cap angular radius (p=3).
    % Higher conflev gives a larger cap (wider angular radius).
    mu = [0 0 1];   kappa = 15;
    for cl = [0.50, 0.90, 0.95, 0.99]
        P = vMFcap(mu, kappa, cl, 1);    % one boundary point is enough
        theta = acosd(mu * P');          % angular radius in degrees
        fprintf('conflev=%.2f: theta_alpha = %.2f deg\n', cl, theta);
    end
%}
%
%{
    % High-dimensional case (p=50): boundary is a Monte Carlo sample.
    % The boundary (p-2)-sphere lives in R^50; view it through 2 coordinates.
    p = 50;
    mu = zeros(1,p); mu(1) = 1;
    kappa = 20;
    P = vMFcap(mu, kappa, 0.95, 300);
    figure; plot(P(:,1),P(:,2),'.');
    xlabel('dim 1'); ylabel('dim 2');
    title(sprintf('p=%d: 2D projection of 300-point cap boundary', p));
%}
%
%{
    % Verify that boundary points have unit norm and lie at angle theta_alpha.
    % For p=3, all points should be exactly on the boundary circle.
    mu = [1 0 0]; kappa = 10; conflev = 0.95;
    P = vMFcap(mu, kappa, conflev, 200);
    norms    = vecnorm(P, 2, 2);           % should all be 1
    cosTheta = P * mu';                    % should all equal t_alpha
    fprintf('norm range: [%.10f, %.10f]\n', min(norms), max(norms));
    fprintf('mu''*P range: [%.10f, %.10f]\n', min(cosTheta), max(cosTheta));
%}

%% Beginning of code

if nargin<4 || isempty(npts)
    npts=100;
end
if nargin<3 || isempty(conflev)
    conflev=0.95;
end

mu=mu(:)';
p=numel(mu);
alphaConf=1-conflev;

if kappa<=0
    talpha=2*alphaConf-1;
elseif p==3
    talpha=1 + log(alphaConf + (1-alphaConf)*exp(-2*kappa))/kappa;
else
    % g is shifted by the constant factor exp(-kappa) (g(t)=exp(kappa*(t-1))*...
    % instead of exp(kappa*t)*...) purely for numerical stability at
    % large kappa: exp(kappa*(t-1))<=1 always for t<=1, so this never
    % overflows, and the constant factor cancels in the normalized CDF.
    g=@(t) exp(kappa*(t-1)) .* (1-t.^2).^((p-3)/2);
    Z=integral(g,-1,1);
    cdf=@(t) integral(g,-1,t)/Z;
    talpha=fzero(@(t) cdf(t)-alphaConf, [-1+1e-10,1-1e-10]);
end
talpha=min(max(talpha,-1),1);
thetaalpha=acos(talpha);

if p==3
    % Exact, deterministic circle
    ref=[1 0 0];
    if abs(mu*ref')>0.9
        ref=[0 1 0];
    end
    e1=cross(ref,mu); e1=e1/norm(e1);
    e2=cross(mu,e1);
    phi=linspace(0,2*pi,npts)';
    P = cos(thetaalpha)*mu + sin(thetaalpha)*(cos(phi)*e1 + sin(phi)*e2);
else
    % Monte Carlo sample of the boundary (p-2)-sphere in the subspace
    % orthogonal to mu: project Gaussian noise off mu, then normalize.
    V=randn(npts,p);
    V=V - (V*mu')*mu;
    normV=sqrt(sum(V.^2,2));
    ok=normV>1e-12;
    V(ok,:)=V(ok,:)./normV(ok);
    if any(~ok)
        % astronomically unlikely (a draw landing exactly along mu);
        % resample those rows
        V(~ok,:)=vMFcap_orthosample(mu,sum(~ok));
    end
    P = cos(thetaalpha)*mu + sin(thetaalpha)*V;
end

end

function V = vMFcap_orthosample(mu,n)
p=numel(mu);
V=zeros(n,p);
for i=1:n
    ok=false;
    while ~ok
        v=randn(1,p);
        v=v-(v*mu')*mu;
        nv=norm(v);
        if nv>1e-12
            V(i,:)=v/nv;
            ok=true;
        end
    end
end
end
%FScategory:CLUS-RobClaMULT
