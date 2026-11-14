function [P, talpha] = vMFcap(mu,kappa,conflev,npts)
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
%   under vMF(mu,kappa) equals conflev, i.e. it returns the boundary of
%
%     {y : mu'*y >= t_alpha}   where P(mu'*y >= t_alpha) = conflev,
%
%   together with the threshold t_alpha.
%
%   For p=3 the boundary is a genuine circle (a 1-sphere embedded in R^3)
%   and vMFcap returns it exactly, as npts evenly spaced points. For p=2
%   the boundary is made of the two end points of an arc. For p>3 the
%   boundary is a (p-2)-sphere in the (p-1)-dimensional subspace
%   orthogonal to mu; vMFcap returns npts points drawn uniformly on that
%   submanifold (in the native R^p) as a Monte Carlo representation. No
%   dimensionality reduction is performed here: to show a cap together
%   with the data use sphereplot (p=2, 3) or sphereprojplot (p>3).
%
%  Required input arguments:
%
%          mu: Mean direction. Vector. 1-by-p unit vector; the centre of
%             the vMF distribution whose confidence cap is requested.
%       kappa: Concentration parameter. Scalar. kappa>=0. kappa=0 gives
%             the uniform distribution on S^(p-1).
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
%              exact circle (the first and the last coincide). For p==2
%              the rows alternate between the two end points of the arc.
%              For p>3 they are a random uniform sample of the boundary
%              (p-2)-sphere in the subspace orthogonal to mu. With npts=0
%              no point is generated (and no random number is used): this
%              is useful when only the threshold talpha is needed.
%                 Example - 200
%                 Data Types - double
%
%  Output:
%
%           P: Boundary points. Matrix. npts-by-p matrix of unit vectors
%              on the boundary of the confidence cap. For p==3 all rows
%              are exactly on the circle in S^2; pass P to sphereplot
%              (option cap). For p>3 P is a Monte Carlo sample of the
%              boundary (p-2)-sphere; see sphereprojplot.
%      talpha: Threshold of the cap. Scalar. The cap is the set of unit
%              vectors y with mu'*y>=talpha, and acos(talpha) is its
%              angular radius.
%
% More About:
%
% ANGULAR RADIUS OF THE CAP
%
% The marginal distribution of t = mu'*y under vMF(mu,kappa) has density
% proportional to (Mardia and Jupp, 2000, Section 9.3):
%
%   g(t) = exp(kappa*t) * (1 - t^2)^((p-3)/2),   t in [-1,1].
%
% The cap threshold t_alpha satisfies P(t >= t_alpha) = conflev, i.e.
% t_alpha is the (1-conflev) quantile of g. The angular radius is then
% theta_alpha = acos(t_alpha).
%
% For p=3 the exponent (p-3)/2 = 0 and g reduces to exp(kappa*t), giving
% the closed form
%
%   t_alpha = 1 + log( (1-conflev) + conflev*exp(-2*kappa) ) / kappa,
%
% which is evaluated with log1p and expm1 to be accurate also for small
% kappa. For kappa=0 and p=3, t is uniform on [-1,1] (Archimedes) and
% t_alpha = 1 - 2*conflev. Note that t is uniform ONLY for p=3: for
% kappa=0 and p~=3 the density of t is (1-t^2)^((p-3)/2) and the general
% procedure is used.
%
% For p=2 the angle theta between y and mu has the von Mises density,
% proportional to exp(kappa*cos(theta)) on (-pi,pi], which is integrated
% numerically; t_alpha is the cosine of the half width of the arc centred
% at mu with probability conflev.
%
% For p>3 there is no closed form and g is integrated numerically. The
% computation is done on the logarithm of g,
%
%   log g(t) = kappa*t + ((p-3)/2)*log(1-t^2),
%
% on a grid which is progressively restricted to the interval where g is
% not negligible. Working with log g is essential in high dimension: for
% p in the hundreds g itself underflows to 0 over the whole interval
% (for example (1-t^2)^((p-3)/2) is of the order of 1e-400 for p=768 and
% t=0.96), and it is concentrated on an interval which is too narrow to
% be found by a general purpose quadrature rule on [-1,1].
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
%   Chichester. (marginal distribution of mu'*y under vMF)
%
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions, Journal of Machine Learning Research, Vol. 6,
%   pp. 1345-1382.
%
% See also: Kentcap, sphereplot, sphereprojplot, vMFrnd, tclustvmf, logvMFpdfFS, kappaFS
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
    %% 95% confidence caps for two vMF distributions on S^2 (p=3).
    % A tight distribution (kappa=40) around [1 0 0] and a diffuse one
    % (kappa=8) around [0 1 0]. The caps are exact circles.
    mu1=[1 0 0]; kappa1=40;
    mu2=[0 1 0]; kappa2=8;
    P1 = vMFcap(mu1,kappa1);
    P2 = vMFcap(mu2,kappa2);
    sphereplot([], 'cap', {P1,P2}, 'center', [mu1; mu2], ...
        'labels', {'kappa=40','kappa=8'});
%}
%
%{
    %% Simulated data with their confidence cap.
    % The cap at level 0.9 is superimposed to 500 units generated from
    % the same distribution; about 90 per cent of them are inside.
    rng(1)
    mu = [0 0 1]; kappa = 15;
    Y = vMFrnd(500, mu, kappa);
    [P, talpha] = vMFcap(mu, kappa, 0.9);
    inside = Y*mu' >= talpha;
    sphereplot(Y, 'group', 2-inside, 'cap', {P}, 'center', mu, ...
        'labels', {'inside','outside'});
    title(sprintf('%.1f%% of the units inside the 90%% cap', 100*mean(inside)))
%}
%
%{
    %% Effect of confidence level on the cap angular radius (p=3).
    % Higher conflev gives a larger cap (wider angular radius).
    mu = [0 0 1];   kappa = 15;
    for cl = [0.50, 0.90, 0.95, 0.99]
        [~, talpha] = vMFcap(mu, kappa, cl, 0);
        fprintf('conflev=%.2f: theta_alpha = %.2f deg\n', cl, acosd(talpha));
    end
%}
%
%{
    %% High-dimensional case (p=300): coverage and projection.
    % The threshold of the cap is computed in dimension 300 and it is
    % checked on simulated data. The data and the cap are then shown with
    % sphereprojplot.
    rng(2)
    p = 300;
    mu = zeros(1,p); mu(1) = 1;
    kappa = 500;
    Y = vMFrnd(5000, mu, kappa);
    [P, talpha] = vMFcap(mu, kappa, 0.95, 300);
    fprintf('Proportion of units inside the 95%% cap: %.4f\n', mean(Y*mu' >= talpha));
    sphereprojplot(Y, 'cap', {P}, 'center', mu);
%}
%
%{
    % Verify that boundary points have unit norm and lie at angle theta_alpha.
    % For p=3, all points should be exactly on the boundary circle.
    mu = [1 0 0]; kappa = 10; conflev = 0.95;
    [P, talpha] = vMFcap(mu, kappa, conflev, 200);
    norms    = vecnorm(P, 2, 2);           % should all be 1
    cosTheta = P * mu';                    % should all equal talpha
    fprintf('norm range: [%.10f, %.10f]\n', min(norms), max(norms));
    fprintf('mu''*P range: [%.10f, %.10f], talpha=%.10f\n', min(cosTheta), max(cosTheta), talpha);
%}
%
%{
    % kappa=0: the uniform distribution in different dimensions.
    % The cap with probability 0.95 covers almost the whole sphere for
    % p=3; in high dimension the uniform distribution concentrates around
    % the "equator" orthogonal to mu and the threshold approaches 0.
    for p = [3 10 100 1000]
        [~, talpha] = vMFcap([1 zeros(1,p-1)], 0, 0.95, 0);
        fprintf('p=%4d: talpha = %.4f\n', p, talpha);
    end
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
kappa=max(kappa,0);

%% Threshold of the cap
if p==3
    if kappa==0
        talpha=1-2*conflev;
    else
        % 1 + log((1-conflev)+conflev*exp(-2*kappa))/kappa
        talpha=1 + log1p(conflev*expm1(-2*kappa))/kappa;
    end
elseif p==2
    % half width a of the arc: the density of the angle from mu is
    % proportional to exp(kappa*cos(theta))
    if kappa==0
        a=pi*conflev;
    else
        thmax=min(pi,12/sqrt(kappa));
        th=linspace(0,thmax,20001);
        F=cumtrapz(th,exp(kappa*(cos(th)-1)));
        F=F/F(end);
        [Fu,iu]=unique(F);
        a=interp1(Fu,th(iu),conflev);
    end
    talpha=cos(a);
else
    % log-density of t=mu'*y, up to a constant
    logg=@(x) kappa*x + ((p-3)/2)*log1p(-x.^2);
    lo=-1; hi=1;
    % the interval is progressively restricted to the region where the
    % density is not negligible
    for it=1:4
        tt=linspace(lo,hi,4001)';
        tt=tt(2:end-1);
        l=logg(tt);
        keep=l>max(l)-60;
        step=tt(2)-tt(1);
        lo=max(lo,min(tt(keep))-step);
        hi=min(hi,max(tt(keep))+step);
    end
    tf=linspace(lo,hi,20001)';
    l=logg(tf);
    l(~isfinite(l))=-Inf;
    F=cumtrapz(tf,exp(l-max(l)));
    F=F/F(end);
    [Fu,iu]=unique(F);
    talpha=interp1(Fu,tf(iu),1-conflev);
end
talpha=min(max(talpha,-1),1);
thetaalpha=acos(talpha);

%% Points on the boundary
if npts==0
    P=zeros(0,p);
elseif p==3
    % Exact, deterministic circle
    ref=[1 0 0];
    if abs(mu*ref')>0.9
        ref=[0 1 0];
    end
    e1=cross(ref,mu); e1=e1/norm(e1);
    e2=cross(mu,e1);
    phi=linspace(0,2*pi,npts)';
    P = cos(thetaalpha)*mu + sin(thetaalpha)*(cos(phi)*e1 + sin(phi)*e2);
elseif p==2
    % the two end points of the arc, alternated
    e1=[-mu(2) mu(1)];
    sgn=(-1).^(0:npts-1)';
    P = cos(thetaalpha)*mu + sin(thetaalpha)*sgn*e1;
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
