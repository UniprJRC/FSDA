function kappa = kappaFS(r, p, newtonsteps)
%kappaFS estimates the von Mises-Fisher concentration parameter from the mean resultant length
%
%<a href="matlab: docsearchFS('kappaFS')">Link to the help page for this function</a>
%
%   kappaFS inverts the equation
%
%       A_p(kappa) = I_{p/2}(kappa) / I_{p/2-1}(kappa) = r
%
%   for kappa>=0, where I_nu is the modified Bessel function of the first
%   kind of order nu, r is the (weighted) mean resultant length of a group
%   of unit vectors in R^p, and p is the ambient dimension (data live on
%   S^(p-1)). A_p is strictly increasing on [0,inf) from 0 to 1, so the
%   equation has a unique solution for every r in [0,1). The function
%   uses the Banerjee et al. (2005) closed-form approximation as the
%   starting value and refines it with Newton-Raphson iterations on A_p
%   (Sra, 2012), bringing the relative error to below 1e-5 uniformly for
%   p in [2,50] and kappa in [0,1000] with just two Newton steps.
%
%   This is the vMF-mixture counterpart of the covariance-matrix update
%   inside the concentration steps of tclust: exactly as tclust recomputes
%   sigmaini(:,:,j) from the current group of untrimmed units, tclustvmf
%   recomputes kappaini(j) by calling kappaFS on the mean resultant length
%   of the units currently assigned to group j.
%
%  Required input arguments:
%
%            r: Mean resultant length(s). Vector or scalar. Values in
%               [0,1); each entry is typically norm(sum(Y_j,1))/n_j for
%               one group of n_j unit vectors (the row-sum of the group
%               matrix normalised by its size). kappaFS is fully
%               vectorized over r, so all groups can be processed at once.
%            p: Dimension of the ambient space. Scalar. Data lie on
%               S^(p-1), so Y has p columns. For the circular case p=2;
%               for the spherical case p=3. p can be large (hundreds or
%               thousands, as for document embeddings).
%  newtonsteps: Number of Newton-Raphson refinement steps. Scalar.
%               Default 2. With newtonsteps=0 kappaFS returns the raw
%               Banerjee et al. (2005) approximation. Each additional
%               step roughly doubles the number of correct significant
%               digits; 2 steps give a relative error below 1e-5
%               uniformly for p in [2,50] and kappa in [0,1000].
%
%  Optional input arguments:
%
% (none beyond what is listed above; all three arguments are positional)
%
%  Output:
%
%        kappa : Estimated concentration parameter(s). Vector or scalar,
%                same size as r. All entries are non-negative and finite:
%                values of r clipped to 1-1e-10 (see More About) are
%                returned as large but finite kappa values.
%
% More About:
%
% APPROXIMATION FORMULA (Banerjee et al., 2005; Sra, 2012)
%
% The starting value is the closed-form approximation:
%
%   kappa_0 = r*(p - r^2) / (1 - r^2)
%
% proposed by Banerjee et al. (2005) and analysed by Sra (2012, eq. 21).
% For r->0 (uniform group) kappa_0->0; for r->1 (degenerate group)
% kappa_0->inf.
%
% NEWTON-RAPHSON REFINEMENT (Sra, 2012)
%
% Each Newton step updates kappa using
%
%   kappa <- kappa - (A_p(kappa) - r) / A_p'(kappa)
%
% where the derivative is (Watson 1944; Sra 2012, eq. 3)
%
%   A_p'(kappa) = 1 - A_p(kappa)^2 - ((p-1)/kappa)*A_p(kappa).
%
% The ratio A_p is computed as exp(log I_{p/2}(kappa)-log I_{p/2-1}(kappa))
% with the two logarithms given by logBesseliFS. In this way the ratio is
% available also when the two Bessel functions overflow (large kappa) or
% underflow (large p and moderate kappa: this is the usual situation for
% embeddings in hundreds of dimensions, where the direct ratio of two
% calls to besseli is 0/0). A safeguard reverts any Newton step that
% produces a non-positive or non-finite result back to the Banerjee
% approximation for that element.
%
% CLIPPING OF r NEAR 1
%
% Values of r very close to 1 (near-degenerate groups) are clipped to
% 1-1e-10 before inversion. This keeps kappa large but finite, which is
% required for the subsequent concentration-parameter restriction step
% restrKappaFS to operate correctly (a kappa of Inf would make the ratio
% max/min kappa undefined and the numerical minimisation inside
% restrKappaFS fail).
%
% References:
%
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions, Journal of Machine Learning Research, Vol. 6,
%   pp. 1345-1382.
%
%   Sra, S. (2012), A short note on parameter approximation for von
%   Mises-Fisher distributions: and a fast implementation of I_s(x),
%   Computational Statistics, Vol. 27, pp. 177-190.
%   doi:10.1007/s00180-011-0232-x
%
%   Watson, G.S. (1983), Statistics on Spheres, Wiley, New York.
%   (derivative of A_p; original result attributed to Watson 1944)
%
% See also: tclustvmf, restrKappaFS, vMFlogC, logBesseliFS, vMFrnd
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('kappaFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Recovery check: for known kappa, verify that kappaFS inverts A_p correctly.
    % For p=5 and several values of kappa, compute r=A_p(kappa) and check
    % that kappaFS(r,p) recovers kappa to within 1e-5 relative error.
    p = 5;
    nu = p/2 - 1;
    Ap = @(k) besseli(nu+1,k,1)./besseli(nu,k,1); % scaled ratio: scaling cancels
    ktrue = [0.5 3 10 50 200];
    r = Ap(ktrue);
    khat = kappaFS(r,p);
    relErr = abs(khat - ktrue)./ktrue;
    disp(array2table([ktrue; khat; relErr],'RowNames',{'kTrue','kHat','relErr'}))
%}
%
%{
    % Effect of newtonsteps: accuracy vs speed trade-off.
    % With newtonsteps=0 the raw Banerjee approximation is returned.
    % Each additional step roughly doubles the significant digits.
    p = 3;
    nu = p/2 - 1;
    Ap = @(k) besseli(nu+1,k,1)./besseli(nu,k,1);
    ktrue = logspace(-1, 3, 50);    % 50 values from 0.1 to 1000
    r = Ap(ktrue);
    for ns = 0:3
        khat = kappaFS(r, p, ns);
        maxRelErr = max(abs(khat - ktrue)./ktrue);
        fprintf('newtonsteps=%d: max relative error = %.2e\n', ns, maxRelErr);
    end
%}
%
%{
    %% Estimation of kappa from simulated vMF data, for several groups at once.
    % kappaFS accepts a vector of mean resultant lengths r (one per group)
    % and returns a vector of kappa estimates in one call.
    rng(1)
    p = 3;
    n1 = 1000; n2 = 1500; n3 = 800;
    ktrue = [40 15 5];
    Y1 = vMFrnd(n1, [1 0 0], ktrue(1));
    Y2 = vMFrnd(n2, [0 1 0], ktrue(2));
    Y3 = vMFrnd(n3, [0 0 1], ktrue(3));
    R = [norm(sum(Y1))/n1, norm(sum(Y2))/n2, norm(sum(Y3))/n3];
    kappa = kappaFS(R, p);
    fprintf('Group mean resultant lengths: %.3f %.3f %.3f\n', R);
    fprintf('True kappa:                   %.2f  %.2f  %.2f\n', ktrue);
    fprintf('Estimated kappa:              %.2f  %.2f  %.2f\n', kappa);
%}
%
%{
    %% High dimension (p=1536, as for document embeddings).
    % For this dimension and this value of kappa the two Bessel functions
    % in A_p(kappa) underflow, and their direct ratio cannot be computed.
    rng(2)
    p = 1536; ktrue = 300;
    Y = vMFrnd(2000, ones(1,p), ktrue);
    r = norm(mean(Y,1));
    fprintf('Direct ratio of Bessel functions at the true kappa: %g\n', ...
        besseli(p/2,ktrue,1)/besseli(p/2-1,ktrue,1));
    fprintf('kappa: true %g, Banerjee approximation %.3f, kappaFS %.3f\n', ...
        ktrue, kappaFS(r,p,0), kappaFS(r,p));
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
    % A_p(kappa) through the logarithms of the two Bessel functions
    Apk = exp(logBesseliFS(nu+1, kappaeval) - logBesseliFS(nu, kappaeval));
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
