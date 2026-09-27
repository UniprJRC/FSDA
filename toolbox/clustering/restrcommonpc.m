function [U, Lambda_vk, info] = restrcommonpc(Sigmaini, niini, ncomp, varargin)
%restrcommonpc forces a common set of principal directions across the k groups
%
%<a href="matlab: docsearchFS('restrcommonpc')">Link to the help function</a>
%
%
%   restrcommonpc implements a "partial common principal components"
%   (PCPC) constraint. It takes the k (unconstrained) group scatter
%   matrices estimated in a tclust M-step and returns, for every group,
%   an orthogonal matrix of eigenvectors U(:,:,j) and a vector of
%   eigenvalues Lambda_vk(:,j) such that the FIRST ncomp columns of
%   U(:,:,j) are IDENTICAL across all k groups (the shared "principal
%   axes"), while the remaining v-ncomp columns/eigenvalues are free and
%   estimated separately for each group.
%
%   Two modes are available, controlled by the optional 'majoraxis'
%   argument (see below):
%
%   * majoraxis=false (the default). The shared axis/axes are chosen as
%     the leading ncomp eigenvector(s) of the (trimming-weighted) pooled
%     scatter matrix M = sum_j (n_j/n) Sigma_j -- a classical, cheap,
%     closed-form, variance-maximizing choice (Ky Fan's theorem, no
%     iteration needed; see "More About" below). This guarantees that
%     the k groups share a common eigenVECTOR, but does NOT guarantee
%     that this shared direction is the DOMINANT (largest-eigenvalue)
%     direction for every group: a group whose own covariance structure
%     disagrees with the pooled vote can end up with the shared axis in
%     one of its minor-axis slots instead.
%
%   * majoraxis=true (only allowed with ncomp=1). The single shared axis
%     is required to be, simultaneously, the MAJOR (dominant,
%     largest-eigenvalue) direction of EVERY group. This is a strictly
%     stronger and generally harder constraint: for v=2 it is solved
%     exactly by a deterministic scan of the constraint circle; for v>2
%     a derivative-free heuristic search on the unit sphere is used
%     instead (see "More About"). A feasible common major axis need not
%     exist (this becomes a real possibility once k>2): when it does
%     not, the direction with the best worst-case margin is returned,
%     and info.feasible is set to false so the caller can tell the two
%     situations apart.
%
%   The output of this function has exactly the same U/Lambda_vk format
%   produced by the usual per-group call
%
%       [Uj,Lambdaj] = eig(Sigmaini(:,:,j));
%       U(:,:,j)     = real(Uj);
%       Lambda_vk(:,j) = real(diag(Lambdaj));
%
%   so it is a drop-in replacement for that loop inside tclust.m: the
%   subsequent eigenvalue restriction (restreigen/restrdeter) and the
%   reconstruction Sigma(:,:,j) = U(:,:,j)*diag(autovalues(:,j))*U(:,:,j)'
%   are completely unaffected and do not need to change, in either mode.
%
%
% Required input arguments:
%
% Sigmaini : unconstrained covariance matrices. v-by-v-by-k array.
%            The k (unconstrained) scatter matrix estimates produced in
%            the current M-step, BEFORE any eigenvalue/shape/determinant
%            restriction is applied.
%               Data Types - single | double
%
%   niini  : size of the groups. Vector of length k.
%            In tclust this already reflects the trimmed subsample (in
%            the crisp case) or the sum of posterior probabilities
%            restricted to the h retained units (mixture case), so
%            trimming is automatically inherited by this routine: units
%            which have been trimmed at the current step play no role in
%            determining the common direction(s).
%               Data Types - single | double
%
%
%  Optional input arguments:
%
%   ncomp  : number of common principal components to impose. Scalar
%            integer, 1 <= ncomp <= v-1. Default is 1 (just the leading
%            axis -- e.g. the major axis of the ellipsoids -- is forced
%            to be the same direction in every group). ncomp=v-1 forces
%            all directions but the last to be common; to force ALL
%            directions to be common (ncomp=v) simply set U(:,:,j) equal
%            for every j using the eigenvectors of the pooled matrix
%            M below -- this provides a cheap, non-iterative, closed
%            form alternative to the "third letter is E" models handled
%            by cpcE/cpcV.
%               Example - 1
%               Data Types - double
%
% The following are given as name-value pairs (varargin):
%
% 'majoraxis' : enforce the shared axis to be each group's own major
%            (dominant) direction. Logical scalar, default false. Only
%            allowed together with ncomp=1 (an error is raised
%            otherwise: "the major axis" only has meaning for a single
%            shared direction). See the description above and "More
%            About" below.
%               Example - 'majoraxis',true
%               Data Types - logical
%
% 'ngrid'  : number of grid points used to scan the circle when v=2 and
%            majoraxis=true. Default 4000 (angular resolution about
%            0.045 degrees). Unused otherwise.
%               Example - 'ngrid',4000
%               Data Types - double
%
% 'niter'  : number of refinement rounds of the general-v heuristic
%            local search used when v>2 and majoraxis=true. Default 200.
%            Unused otherwise (v=2 uses the exact grid search; when
%            majoraxis=false no search is needed at all).
%               Example - 'niter',200
%               Data Types - double
%
%
% Output:
%
%        U : eigenvectors. v-by-v-by-k array.
%            U(:,:,j) is an orthogonal matrix. Its first ncomp columns
%            are identical for every j=1,...,k (the common principal
%            axes); the columns from ncomp+1 to v are group specific.
%
% Lambda_vk : eigenvalues. v-by-k matrix.
%            Column j contains the eigenvalues associated with the
%            columns of U(:,:,j), in the same order. These are NOT yet
%            restricted: they still have to go through
%            restreigen/restrdeter exactly as in the unconstrained case.
%
%     info : diagnostic structure with fields:
%            info.feasible : logical.
%                 majoraxis=false: always true (a shared eigenvector
%                    always exists by construction).
%                 majoraxis=true: true if a direction was found that is
%                    simultaneously the dominant direction of every
%                    group; false if the returned axis is only the best
%                    achievable worst-case compromise (see "More About").
%            info.margin   : k-by-1 vector, computed whenever ncomp=1
%                 (empty otherwise). margin(j) = (variance of group j
%                 along the returned shared axis) minus (its largest
%                 variance among directions orthogonal to the shared
%                 axis). margin(j)>=0 means the shared axis really is
%                 group j's own major axis; margin(j)<0 means a minor
%                 axis was forced into the "major" slot for that group.
%                 This diagnostic is computed and returned in BOTH
%                 modes, so it can be used to check, after the fact,
%                 whether the default (majoraxis=false) mode happened to
%                 produce a shared axis that is every group's major axis
%                 or not.
%            info.u        : the returned shared unit vector, v-by-1
%                 (only when ncomp=1; empty otherwise).
%
%
% More About:
%
% Given weights $n_j$ (already reflecting trimming) and unconstrained
% scatter matrices $\Sigma_j$, $j=1,\ldots,k$, define the pooled matrix
%
% \[
% M = \sum_{j=1}^k \frac{n_j}{n} \Sigma_j , \qquad n=\sum_j n_j .
% \]
%
% DEFAULT MODE (majoraxis=false). For a single common direction, the
% unit vector $u$ which maximizes the pooled (trimmed) explained
% variance $\Phi(u)=\sum_j n_j\, u'\Sigma_j u = u' M u$ is simply the
% eigenvector of $M$ associated with its largest eigenvalue: $\Phi$ is
% linear in the $\Sigma_j$'s, so the maximization reduces to an ordinary
% Rayleigh-quotient problem for $M$. By the same (Ky Fan) argument, the
% best set of ncomp mutually orthogonal common directions is simply
% given by the leading ncomp eigenvectors of $M$ -- no iterative
% algorithm is required. This criterion maximizes TOTAL pooled variance
% explained, and can therefore be dominated by whichever group has the
% largest scale/sample size, at the expense of a poor fit for the
% others; it does not by itself guarantee that the shared axis is any
% individual group's own major axis.
%
% MAJORAXIS MODE (majoraxis=true, ncomp=1). Here a direction $u$ is
% sought such that $u'\Sigma_j u \ge w'\Sigma_j w$ for every unit vector
% $w\perp u$ and every group $j$ -- i.e. $u$ is simultaneously the
% dominant eigenvector of every $\Sigma_j$. For $v=2$, writing each
% group's own eigen-decomposition as eigenvalues
% $\lambda_{j1}\ge\lambda_{j2}$ at angle $\theta_j$, the set of angles
% $\phi$ for which $\phi$ is group $j$'s OWN major axis is exactly the
% 90-degree-wide arc $|\phi-\theta_j|<45^\circ$ (mod $180^\circ$).
% Because any two undirected axes in the plane are at most $90^\circ$
% apart, for k=2 groups a feasible common major axis exists for almost
% any configuration (the only knife-edge case is two groups whose own
% axes are exactly $90^\circ$ apart, where the two arcs meet only at
% their shared boundary with zero margin). For k>2 groups whose own
% orientations span more than $90^\circ$, no common major axis can
% exist. This function scans $\phi\in[0,\pi)$ on a fine grid; among
% angles for which every group's margin is non-negative (feasible), it
% returns the one maximizing the pooled-variance criterion above;
% otherwise it returns the angle maximizing the worst-case (minimum
% over groups) margin, and sets info.feasible=false. For v>2 there is no
% such closed-form scan; a heuristic derivative-free search on the unit
% sphere is used instead, seeded from each group's own leading
% eigenvector and from the default-mode (pooled) direction (see the
% local function majoraxisND below); it is not guaranteed to find the
% global optimum for v>2, unlike the v=2 case which is exact.
%
%
% See also: tclust, restreigen, restrdeter, restrSigmaGPCM
%
%
% References:
%
%   Flury B.D. (1988), "Common Principal Components and Related
%   Multivariate Models", Wiley.
%
%   Garcia-Escudero L.A., Mayo-Iscar, A. and Riani M. (2020). Model-based
%   clustering with determinant-and-shape constraint, Statistics and
%   Computing, vol. 30, pp. 1363-1380.
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('restrcommonpc')">Link to the help function</a>
%
%$LastChangedDate:: $: Date of the last commit

% Examples:

%{
    % Two groups, p=3, force the major axis to be common (default mode).
    rng(1)
    v=3; k=2;
    A1=randn(v); Sigma1=A1*A1';
    A2=randn(v); Sigma2=A2*A2';
    Sigmaini=cat(3,Sigma1,Sigma2);
    niini=[120;80];
    [U,Lambda_vk]=restrcommonpc(Sigmaini,niini,1);
    % The first column of U(:,:,1) and U(:,:,2) must coincide (up to sign)
    disp(abs(U(:,1,1))-abs(U(:,1,2)))
%}

%{
    % Same data, but require the shared axis to be each group's own
    % major axis (majoraxis mode), and inspect the feasibility diagnostic.
    rng(1)
    v=3; k=2;
    A1=randn(v); Sigma1=A1*A1';
    A2=randn(v); Sigma2=A2*A2';
    Sigmaini=cat(3,Sigma1,Sigma2);
    niini=[120;80];
    [U,Lambda_vk,info]=restrcommonpc(Sigmaini,niini,1,'majoraxis',true);
    disp(info.feasible)
    disp(info.margin)
%}

%% Beginning of code
[v,~,k] = size(Sigmaini);

if nargin<3 || isempty(ncomp)
    ncomp = 1;
end
ncomp = round(ncomp);
if ncomp<1
    ncomp = 1;
end
if ncomp>v-1
    ncomp = v-1;
end

options = struct('majoraxis',false,'ngrid',4000,'niter',200);
if nargin>3
    UserOptions = varargin(1:2:length(varargin));
    if ~isempty(UserOptions)
        checkoptions = cellfun(@(x) any(strcmpi(x,fieldnames(options))), UserOptions);
        if any(~checkoptions)
            error('FSDA:restrcommonpc:WrongInputOpt','Unrecognized option in restrcommonpc');
        end
        for i=1:2:length(varargin)
            options.(varargin{i}) = varargin{i+1};
        end
    end
end
majoraxis = options.majoraxis;
ngrid = options.ngrid;
niter = options.niter;

if majoraxis && ncomp~=1
    error('FSDA:restrcommonpc:WrongInputOpt', ...
        '''majoraxis'',true requires ncomp=1: "the major axis" only has meaning for a single shared direction.');
end

niini = niini(:);
sumnini = sum(niini);
w = niini/sumnini;

% Pooled, trimming-weighted scatter matrix. Groups with niini(j)==0
% (e.g. temporarily empty clusters) simply do not contribute.
M = zeros(v,v);
for j=1:k
    if niini(j) > 0
        M = M + w(j) * Sigmaini(:,:,j);
    end
end
% enforce exact symmetry to avoid spurious complex eigenvalues from
% rounding
M = (M+M')/2;

if majoraxis
    % ---- majoraxis mode: shared direction must be each group's own major axis
    if v==2
        [u, margin_j, feasible] = majoraxis2D(Sigmaini, w, ngrid);
    else
        [u, margin_j, feasible] = majoraxisND(Sigmaini, w, niter);
    end
    Ucommon = u;
    B = null(u');
else
    % ---- default mode: leading ncomp eigenvectors of the pooled matrix
    [Vall, Dall] = eig(M);
    dM = real(diag(Dall));
    [~, ord] = sort(dM, 'descend');
    Vall = real(Vall(:,ord));

    Ucommon = Vall(:, 1:ncomp);        % shared directions (same for all j)
    B       = Vall(:, ncomp+1:end);    % fixed orthonormal basis of their
                                        % orthogonal complement (v-by-(v-ncomp))
    feasible = true;
    margin_j = [];
    if ncomp==1
        % diagnostic only: is the (single) chosen shared axis actually
        % each group's own major axis, even though this was not enforced?
        margin_j = zeros(k,1);
        for j=1:k
            Sj = Sigmaini(:,:,j); Sj=(Sj+Sj')/2;
            a = Ucommon'*Sj*Ucommon;
            if v>1
                Rj = B'*Sj*B; Rj=(Rj+Rj')/2;
                lam2max = max(eig(Rj));
            else
                lam2max = -Inf;
            end
            margin_j(j) = a-lam2max; %#ok<AGROW>
        end
        feasible = all(margin_j >= -1e-8);
    end
end

U = NaN(v,v,k);
Lambda_vk = NaN(v,k);

for j=1:k
    Sj = Sigmaini(:,:,j);
    Sj = (Sj+Sj')/2;

    % Variance of group j along each of the ncomp shared axes
    lamShared = diag(Ucommon' * Sj * Ucommon);

    if ncomp < v
        Rj = B' * Sj * B;
        Rj = (Rj+Rj')/2;
        [Wj, Dj] = eig(Rj);
        Wj = real(Wj);
        dj = real(diag(Dj));
        [dj, ord2] = sort(dj,'descend');
        Wj = Wj(:,ord2);
        Uresidual = B*Wj;   % v-by-(v-ncomp), orthogonal to Ucommon
    else
        dj = zeros(0,1);
        Uresidual = zeros(v,0);
    end

    U(:,:,j) = [Ucommon, Uresidual];
    Lambda_vk(:,j) = [lamShared; dj];
end

info.feasible = feasible;
info.margin = margin_j;
if ncomp==1
    info.u = Ucommon;
else
    info.u = [];
end

end

%% ------------------------------------------------------------------
function [u, margin_j, feasible] = majoraxis2D(Sigmaini, w, ngrid)
% Exact, deterministic solution for v=2 via a fine grid scan over the
% (pi-periodic) angle of the shared axis, restricted to the feasible
% region whenever one exists (see "More About" in the main help).
k = size(Sigmaini,3);
phis = linspace(0,pi,ngrid+1);
phis(end) = [];

marginmat = zeros(numel(phis),k);
objvec = zeros(numel(phis),1);
for idx = 1:numel(phis)
    phi = phis(idx);
    uu = [cos(phi); sin(phi)];
    ww = [-sin(phi); cos(phi)];
    s = 0;
    for j=1:k
        Sj = Sigmaini(:,:,j);
        a = uu'*Sj*uu;
        b = ww'*Sj*ww;
        marginmat(idx,j) = a-b;
        s = s + w(j)*a;
    end
    objvec(idx) = s;
end

feasmask = all(marginmat>=0,2);
if any(feasmask)
    idxfeas = find(feasmask);
    [~,best] = max(objvec(idxfeas));
    bestidx = idxfeas(best);
    feasible = true;
else
    worst = min(marginmat,[],2);
    [~,bestidx] = max(worst);
    feasible = false;
end

phistar = phis(bestidx);
u = [cos(phistar); sin(phistar)];
margin_j = marginmat(bestidx,:)';
end

%% ------------------------------------------------------------------
function [u, margin_j, feasible] = majoraxisND(Sigmaini, w, niter)
% Heuristic solution for v>2: seed candidates with each group's own
% leading eigenvector plus the pooled Ky-Fan direction, keep the
% candidate with the best worst-case margin, then refine with a simple
% derivative-free (shrinking-step, random-tangent) local search on the
% sphere. Not guaranteed globally optimal; for a rigorous general-v
% solution consider a proper manifold-optimization toolbox.
[v,~,k] = size(Sigmaini);

worstmargin = @(uu) minmargin(uu,Sigmaini);

cands = zeros(v,k+1);
for j=1:k
    [Vj,Dj] = eig(Sigmaini(:,:,j));
    [~,imax] = max(diag(Dj));
    cands(:,j) = Vj(:,imax);
end
M = zeros(v,v);
for j=1:k
    M = M + w(j)*Sigmaini(:,:,j);
end
M = (M+M')/2;
[VM,DM] = eig(M);
[~,imax] = max(diag(DM));
cands(:,k+1) = VM(:,imax);

bestval = -Inf;
bestu = cands(:,1);
for c = 1:size(cands,2)
    val = worstmargin(cands(:,c));
    if val>bestval
        bestval = val;
        bestu = cands(:,c);
    end
end

u = bestu/norm(bestu);
step = 0.2;
for it = 1:niter
    improved = false;
    for trial = 1:2*v
        d = randn(v,1);
        d = d - (d'*u)*u;
        if norm(d)<eps
            continue
        end
        d = d/norm(d);
        utry = u + step*d;
        utry = utry/norm(utry);
        val = worstmargin(utry);
        if val > bestval + 1e-10
            bestval = val;
            u = utry;
            improved = true;
        end
    end
    if ~improved
        step = step/2;
        if step<1e-6
            break
        end
    end
end

margin_j = zeros(k,1);
for j=1:k
    Sj = Sigmaini(:,:,j);
    a = u'*Sj*u;
    B = null(u');
    R = B'*Sj*B;
    lam2max = max(eig((R+R')/2));
    margin_j(j) = a-lam2max;
end
feasible = all(margin_j >= -1e-8);
end

function m = minmargin(u, Sigmaini)
u = u/norm(u);
k = size(Sigmaini,3);
B = null(u');
vals = zeros(k,1);
for j=1:k
    Sj = Sigmaini(:,:,j);
    a = u'*Sj*u;
    R = B'*Sj*B;
    lam2max = max(eig((R+R')/2));
    vals(j) = a-lam2max;
end
m = min(vals);
end
%FScategory:CLUS-RobClaMULT
