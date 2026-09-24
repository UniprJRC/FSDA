function [U, Lambda_vk] = restrcommonpc(Sigmaini, niini, ncomp)
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
%   The output of this function has exactly the same format produced by
%   the usual per-group call
%
%       [Uj,Lambdaj] = eig(Sigmaini(:,:,j));
%       U(:,:,j)     = real(Uj);
%       Lambda_vk(:,j) = real(diag(Lambdaj));
%
%   so it is a drop-in replacement for that loop inside tclust.m: the
%   subsequent eigenvalue restriction (restreigen/restrdeter) and the
%   reconstruction Sigma(:,:,j) = U(:,:,j)*diag(autovalues(:,j))*U(:,:,j)'
%   are completely unaffected and do not need to change.
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
% For a single common direction, the unit vector $u$ which maximizes the
% pooled (trimmed) explained variance
% $\Phi(u)=\sum_j n_j\, u'\Sigma_j u = u' M u$
% is simply the eigenvector of $M$ associated with its largest
% eigenvalue: $\Phi$ is linear in the $\Sigma_j$'s, so the maximization
% reduces to an ordinary Rayleigh-quotient problem for $M$. By the same
% (Ky Fan) argument, the best set of ncomp mutually orthogonal common
% directions is simply given by the leading ncomp eigenvectors of $M$ --
% no iterative algorithm is required.
%
% Once the common axes $U_{common}$ are fixed, for every group $j$ we:
%   (a) keep $a_{j} = \mathrm{diag}(U_{common}' \Sigma_j U_{common})$ as
%       the (group specific) variances along the shared directions;
%   (b) project $\Sigma_j$ onto the (fixed) orthogonal complement $B$ of
%       $U_{common}$, i.e. $R_j = B' \Sigma_j B$, and diagonalize $R_j$
%       to get the group-specific remaining eigenvectors/eigenvalues.
% This is equivalent to replacing $\Sigma_j$ by its closest (in the
% sense of dropping the cross terms between the shared subspace and its
% complement) covariance matrix which has $U_{common}$ as a set of exact
% eigenvectors. The construction guarantees each returned matrix
% $U(:,:,j)\,diag(Lambda\_vk(:,j))\,U(:,:,j)'$ is symmetric PSD.
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
    % Two groups, p=3, force the major axis to be common.
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

niini = niini(:);
sumnini = sum(niini);

% Pooled, trimming-weighted scatter matrix. Groups with niini(j)==0
% (e.g. temporarily empty clusters) simply do not contribute.
M = zeros(v,v);
for j=1:k
    if niini(j) > 0
        M = M + (niini(j)/sumnini) * Sigmaini(:,:,j);
    end
end
% enforce exact symmetry to avoid spurious complex eigenvalues from
% rounding
M = (M+M')/2;

[Vall, Dall] = eig(M);
dM = real(diag(Dall));
[~, ord] = sort(dM, 'descend');
Vall = real(Vall(:,ord));

Ucommon = Vall(:, 1:ncomp);        % shared directions (same for all j)
B       = Vall(:, ncomp+1:end);    % fixed orthonormal basis of their
                                    % orthogonal complement (v-by-(v-ncomp))

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
        Uresidual = B*Wj;   % v-by-(v-ncomp), orthogonal to Ucommon
    else
        dj = zeros(0,1);
        Uresidual = zeros(v,0);
    end

    U(:,:,j) = [Ucommon, Uresidual];
    Lambda_vk(:,j) = [lamShared; dj];
end

end
%FScategory:CLUS-RobClaMULT
