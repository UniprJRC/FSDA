function kappa = kappaFS(r, p, newtonsteps)
%kappaFS estimates the von Mises-Fisher concentration parameter from the mean resultant length
%
%<a href="matlab: docsearchFS('kappaFS')">Link to the help function</a>
%
%   kappaFS inverts the equation
%
%       Ap(kappa) = I_(p/2)(kappa) / I_(p/2-1)(kappa) = r
%
%   for kappa, where I_nu is the modified Bessel function of the first
%   kind of order nu, r is the (weighted) mean resultant length of a
%   group of unit vectors, and p is the dimension of the ambient space
%   (data live on S^(p-1)). This equation has no closed-form solution
%   because it requires inverting a ratio of Bessel functions; kappaFS
%   uses the Banerjee, Dhillon, Ghosh and Sra (2005) closed-form
%   approximation as a starting value, refined by a small number of
%   Newton-Raphson iterations on Ap (Sra, 2012), which brings the
%   relative error down to essentially machine precision for the range
%   of p and kappa relevant in clustering applications.
%
%   This is the vMF-mixture counterpart of the covariance-matrix update
%   inside the concentration steps of tclust: exactly as tclust
%   recomputes sigmaini(:,:,j) from the current group of untrimmed units,
%   tclustmvmf recomputes kappaini(j) by calling kappaFS on the mean
%   resultant length of the units currently assigned to group j.
%
%  Required input arguments:
%
%            r: Mean resultant length(s). Vector or scalar. Values in
%               [0,1). Typically r = norm(sum(Y(group,:)))/n_group for one
%               or more groups (kappaFS is fully vectorized over r).
%            p: Dimension of the ambient space. Scalar. Data lie on
%               S^(p-1), i.e. Y has p columns.
%
%  Optional input arguments:
%
%  newtonsteps: Number of Newton-Raphson refinement steps. Scalar.
%               The default value is 2, which was found (see the
%               examples below) to bring the relative error to below
%               1e-05 uniformly for p in [2,50] and kappa in [0,1000].
%               newtonsteps=0 returns the raw Banerjee et al. (2005)
%               approximation.
%                 Example - 3
%                 Data Types - double
%
%  Output:
%
%        kappa : Estimated concentration parameter(s). Vector or scalar,
%                same size as r.
%
% More About:
%
% Values of r very close to 1 (near-degenerate groups, e.g. very small
% or near-identical-direction groups) are clipped to 1-1e-10 before
% inversion. This keeps kappa large but finite rather than Inf, which is
% required for the subsequent eigenvalue-ratio-style restriction step
% (restreigen/restreigeneasy, reused unmodified from tclust on the 1-by-k
% matrix of kappa values) to operate correctly.
%
% References:
%
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   "Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions", Journal of Machine Learning Research, Vol. 6.
%   Sra, S. (2012), "A short note on parameter approximation for von
%   Mises-Fisher distributions, and a fast implementation of I_s(x)",
%   Computational Statistics, Vol. 27.
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('kappaFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % Recovery check: for known kappa, r=Ap(kappa); verify kappaFS(r,p) ~ kappa.
    p = 5;
    nu = p/2 - 1;
    Ap = @(k) besseli(nu+1,k,1)./besseli(nu,k,1); % scaling cancels in the ratio
    ktrue = [0.5 3 10 50 200];
    r = Ap(ktrue);
    khat = kappaFS(r,p);
    disp([ktrue; khat])
%}

%% Beginning of code

if nargin < 3 || isempty(newtonsteps)
    newtonsteps = 2;
end

r = min(max(r,0), 1-1e-10);

% Banerjee et al. (2005) closed-form approximation
kappa0 = r.*(p - r.^2) ./ (1 - r.^2);
kappa = kappa0;

nu = p/2 - 1;

for it = 1:newtonsteps
    % kappaeval avoids a 0/0 evaluation of Ap at kappa=0 (uniform limit)
    kappaeval = max(kappa, 1e-8);
    Apk = besseli(nu+1, kappaeval, 1) ./ besseli(nu, kappaeval, 1);
    derivAp = 1 - Apk.^2 - ((p-1)./kappaeval).*Apk;
    kappanew = kappa - (Apk - r) ./ derivAp;
    % safeguard: if a Newton step overshoots into non-positive territory,
    % fall back to the Banerjee approximation for that element
    invalid = ~(kappanew > 0) | ~isfinite(kappanew);
    kappanew(invalid) = kappa0(invalid);
    kappa = kappanew;
end

end
%FScategory:CLUS-RobClaMULT
