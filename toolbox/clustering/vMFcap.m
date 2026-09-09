function P = vMFcap(mu,kappa,conflev,npts)
%vMFcap computes the boundary of a von Mises-Fisher confidence cap in any dimension
%
%<a href="matlab: docsearchFS('vMFcap')">Link to the help function</a>
%
%   vMFcap is the directional-data counterpart of a Gaussian confidence
%   ellipse: because the vMF(mu,kappa) density depends on y only through
%   the angle between y and mu, its iso-density contours are exactly the
%   spherical caps {y : angle(y,mu) <= theta}, for any embedding
%   dimension p (data on S^(p-1) in R^p). vMFcap returns npts points on
%   the boundary of the cap {angle(y,mu) <= theta_alpha}, theta_alpha
%   chosen so the cap has probability conflev under vMF(mu,kappa).
%
%   For p=3 the boundary is a genuine circle (a 1-sphere) and vMFcap
%   returns it exactly, deterministically parametrized. For p~=3 the
%   boundary is a (p-2)-sphere living in the (p-1)-dimensional subspace
%   orthogonal to mu -- a whole submanifold, not a curve -- so vMFcap
%   instead returns npts points drawn (uniformly on that submanifold) as
%   a Monte Carlo representation of the boundary, still living in the
%   native R^p (no dimensionality reduction is done here: choosing how
%   to view a high-dimensional cap in 2 or 3 dimensions is left to the
%   caller, e.g. tclustmvmf's 'plotdims' option).
%
%  Required input arguments:
%
%          mu: Mean direction. Vector. 1-by-p unit vector.
%       kappa: Concentration parameter. Scalar. kappa>=0.
%
%  Optional input arguments:
%
%     conflev: Confidence level of the cap. Scalar in (0,1). Default 0.95.
%                 Example - 0.90
%       npts : Number of points used to represent the boundary. Scalar.
%              Default 100. For p==3 these are evenly spaced around the
%              exact circle; for p~=3 they are a random sample of the
%              boundary (p-2)-sphere.
%                 Example - 200
%
%  Output:
%
%           P: npts-by-p matrix of unit vectors on the boundary of the
%              confidence cap. For p==3, suitable directly for
%              plot3(P(:,1),P(:,2),P(:,3)). For p~=3, project P (and the
%              data) onto whatever 2D/3D view is of interest before
%              plotting.
%
% More About:
%
% The angular radius theta_alpha is found from the marginal distribution
% of t=cos(angle(y,mu)) under vMF(mu,kappa), which (regardless of p) has
% density proportional to
%
%       g(t) = exp(kappa*t) * (1-t^2)^((p-3)/2),   t in [-1,1]
%
% (a consequence of the surface measure on S^(p-1) factoring into a
% uniform azimuthal part and this "polar angle" part). For p=3 the
% exponent is 0 and g reduces to a simple exponential tilt, giving the
% closed form
%
%       t_alpha = 1 + log( (1-conflev) + conflev*exp(-2*kappa) ) / kappa
%
% used here whenever p==3 (rewritten to avoid overflow in sinh(kappa) at
% large kappa). For p~=3 there is no closed form (g involves a
% Beta-type, not exponential-family-linear, term), so t_alpha is instead
% found by directly integrating g numerically and inverting the
% resulting CDF with a root-finder. For kappa=0 (uniform on the sphere,
% any p) t_alpha=2*conflev-1 exactly, since t is then uniform on [-1,1]
% regardless of p.
%
% Once theta_alpha=acos(t_alpha) is known, any point on the boundary is
% cos(theta_alpha)*mu + sin(theta_alpha)*v for v a unit vector in the
% (p-1)-dimensional subspace orthogonal to mu. For p==3 that subspace is
% 2-dimensional and v=[cos(phi) sin(phi)] in an orthonormal basis of it
% traces the exact circle. For general p, v is instead drawn at random
% (uniformly) from the (p-2)-sphere of unit vectors in that subspace.
%
% See also: logvMFpdfFS, kappaFS, tclustmvmf
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('vMFcap')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % 95% confidence cap around a tight and a diffuse mean direction, p=3.
    mu1=[1 0 0]; kappa1=40;
    mu2=[0 1 0]; kappa2=8;
    P1 = vMFcap(mu1,kappa1);
    P2 = vMFcap(mu2,kappa2);
    [Yx,Yy,Yz]=sphere(40);
    figure; surf(Yx,Yy,Yz,'FaceColor',[0.9 0.9 1],'EdgeColor','none','FaceAlpha',0.25);
    hold on
    plot3(P1(:,1),P1(:,2),P1(:,3),'b-','LineWidth',2)
    plot3(P2(:,1),P2(:,2),P2(:,3),'r-','LineWidth',2)
    axis equal; view(3)
%}

%{
    % High-dimensional case (p=50, e.g. embeddings): the boundary is a
    % Monte Carlo sample of a 48-sphere living in R^50; view it through
    % an arbitrary pair of coordinates.
    p=50;
    mu=zeros(1,p); mu(1)=1;
    kappa=20;
    P = vMFcap(mu,kappa,0.95,300);
    plot(P(:,1),P(:,2),'.')
    xlabel('dim 1'); ylabel('dim 2')
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
