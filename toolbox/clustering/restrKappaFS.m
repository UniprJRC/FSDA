function kappa = restrKappaFS(kappaHat, r, niini, p, restrfactor, tol)
%restrKappaFS restricts the vMF concentration parameters across groups
%
%<a href="matlab: docsearchFS('restrKappaFS')">Link to the help function</a>
%
%   restrKappaFS is the von Mises-Fisher analogue of restreigen: given
%   raw (unconstrained MLE) concentration parameters kappaHat_j for each
%   of k groups, it returns restricted values with
%   max(kappa)/min(kappa) <= restrfactor, chosen to maximize the total
%   (niini-weighted) vMF log-likelihood subject to that constraint.
%
%   Like restreigen, this follows Fritz, Garcia-Escudero and
%   Mayo-Iscar (2013): for a fixed "floor" m, the constrained per-group
%   optimum is the raw estimate clipped to [m, restrfactor*m], because the
%   per-group loss
%
%       L_j(kappa) = -log Cp(kappa) - kappa*r_j
%
%   (minus the average vMF log-likelihood of group j, r_j being its mean
%   resultant length) is convex in kappa -- its second derivative equals
%   Ap'(kappa) > 0, where Ap is the Bessel-ratio function inverted by
%   kappaFS -- with unconstrained minimizer exactly kappaHat_j =
%   kappaFS(r_j,p). Unlike the Gaussian eigenvalue case handled by
%   restreigen, where the analogous per-direction loss log(e)+d/e is
%   simple enough that the optimal shared floor m has a closed-form
%   algebraic solution, log Cp(kappa) involves a modified Bessel
%   function of kappa, so there is no closed form here: restrKappaFS
%   finds the optimal m by direct 1-D numerical minimization of the
%   total weighted loss over log(m) (log-space keeps the search
%   well-scaled even when a near-singleton group pushes the raw kappaHat
%   range across many orders of magnitude), rather than by the
%   breakpoint algebra used inside restreigen.
%
%  Required input arguments:
%
%    kappaHat: Raw (unrestricted) concentration MLEs. Vector. k-by-1,
%              typically kappaFS(r,p).
%           r: Mean resultant lengths. Vector. k-by-1, one per group,
%              used (together with p) to evaluate the vMF loss; NOT
%              recoverable from kappaHat alone without redundant Bessel
%              inversion, hence passed separately.
%       niini: Cluster sizes. Vector. k-by-1. Groups with niini=0 (or
%              kappaHat=NaN) are treated as empty: they are excluded from
%              the restriction and, on return, set to the mean of the
%              restricted values of the non-empty groups (mirroring
%              restreigen's own convention for empty groups).
%           p: Dimension of the ambient space. Scalar.
% restrfactor: Restriction factor. Scalar >=1.
%
%  Optional input arguments:
%
%         tol: Tolerance passed to fminbnd. Scalar. Default 1e-08.
%
%  Output:
%
%       kappa: Restricted concentration parameters. Vector, k-by-1, same
%              size as kappaHat. max(kappa)/min(kappa) <= restrfactor
%              among the originally non-empty groups.
%
% See also: restreigen, kappaFS, vMFlogC, tclustmvmf
%
% References:
%
%   Fritz H., Garcia-Escudero, L.A. and Mayo-Iscar, A. (2013), A fast
%   algorithm for robust constrained clustering, Computational
%   Statistics and Data Analysis, Vol. 61, pp. 124-136.
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('restrKappaFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    % A tight, spurious near-singleton group (kappa artificially huge)
    % gets pulled back towards the other groups.
    kappaHat = [8000; 25; 30];
    r        = [0.9999; 0.85; 0.87];
    niini    = [3; 140; 150];
    p = 3;
    kappa = restrKappaFS(kappaHat,r,niini,p,20);
    disp(kappa')
    disp(max(kappa)/min(kappa))
%}

%% Beginning of code

if nargin<6 || isempty(tol)
    tol=1e-08;
end

kappa = kappaHat;

ok = ~isnan(kappaHat) & niini>0;
if sum(ok) <= 1
    if any(~ok) && any(ok)
        kappa(~ok) = kappa(ok);
    end
    return
end

kok = kappaHat(ok);
rok = r(ok);
nok = niini(ok);
w   = nok / sum(nok);

if max(kok)/min(kok) <= restrfactor
    if any(~ok)
        kappa(~ok) = mean(kok);
    end
    return
end

c = restrfactor;

    function s = totalloss(logm)
        m = exp(logm);
        e = min(max(kok, m), c*m);
        s = sum( w .* ( -vMFlogC(e,p) - e.*rok ) );
    end

logmlow  = log(min(kok)/c);
logmhigh = log(max(kok));
logmopt = fminbnd(@totalloss, logmlow, logmhigh, optimset('TolX',tol,'Display','off'));
mopt = exp(logmopt);

krestr = min(max(kok, mopt), c*mopt);
kappa(ok) = krestr;
if any(~ok)
    kappa(~ok) = mean(krestr);
end

end
%FScategory:CLUS-RobClaMULT
