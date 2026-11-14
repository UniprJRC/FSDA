function [kappaout, betaout] = restrKentFS(kappa, beta, r, dq, ni, restrfactor, tol)
%restrKentFS applies the concentration restriction to a set of Kent mixture components
%
%<a href="matlab: docsearchFS('restrKentFS')">Link to the help page for this function</a>
%
%   restrKentFS is the Kent-distribution analogue of restrKappaFS: it
%   enforces the constraint
%
%       max_j(kappa_j) / min_j(kappa_j)  <=  restrfactor
%
%   across the groups, in the way which is optimal for the Kent
%   likelihood. When the concentration kappa_j of a group has to be
%   changed, its ovalness beta_j is rescaled proportionally, so that the
%   ratio beta_j/kappa_j (which determines the shape, i.e. the
%   eccentricity, of the contours of the group) is preserved and the
%   condition beta_j < kappa_j/2 keeps holding. The restriction therefore
%   acts on the "size" of the groups and leaves their shape and
%   orientation untouched.
%
%   If the unrestricted estimates already satisfy the constraint they are
%   returned unchanged. Otherwise, exactly as in restrKappaFS, the
%   restricted estimates are the unrestricted ones truncated to an
%   interval [m, restrfactor*m], and the optimal m is found by a 1-D
%   minimisation of the (minus) log-likelihood (see "More About").
%
%   When restrfactor=1 all kappa_j are set to the same value. Groups with
%   NaN kappa (empty/degenerate groups, as flagged by KentMstepFS) are
%   passed through unchanged and do not take part in the restriction.
%
%  Required input arguments:
%
%       kappa: Concentration parameters. k-by-1 vector. Unrestricted
%              estimates from KentMstepFS. Entries may be NaN (empty
%              group).
%        beta: Ovalness parameters. k-by-1 vector. Unrestricted estimates
%              from KentMstepFS.
%           r: Mean resultant lengths. k-by-1 vector. r(j) is the
%              (weighted) mean of gamma3_j'*y over the units of group j,
%              as returned by KentMstepFS (6th output).
%          dq: Differences of second moments. k-by-1 vector. dq(j) is the
%              (weighted) mean of (gamma1_j'*y)^2-(gamma2_j'*y)^2 over
%              the units of group j, as returned by KentMstepFS (7th
%              output). r and dq are the sufficient statistics of
%              (kappa_j,beta_j) for given orientation.
%          ni: Group sizes. k-by-1 vector. Number of units (or sum of the
%              posterior probabilities) of each group. They are the
%              weights of the groups in the likelihood.
% restrfactor: Restriction factor. Scalar >= 1. Maximum allowed ratio
%              max(kappa_j)/min(kappa_j). restrfactor=1 forces a common
%              kappa.
%
%  Optional input arguments:
%
%          tol: Relative tolerance for the 1-D minimisation. Scalar.
%               Default 1e-8.
%                 Example - 1e-10
%                 Data Types - double
%
%  Output:
%
%     kappaout: Restricted concentration parameters. k-by-1 vector, with
%               max(kappaout)/min(kappaout) <= restrfactor. NaN for
%               empty/degenerate groups.
%      betaout: Restricted ovalness parameters. k-by-1 vector.
%               betaout(j) = kappaout(j)*beta(j)/kappa(j).
%
% More About:
%
% LOSS FUNCTION
%
% For given orientation, minus the log-likelihood of group j divided by
% n_j is
%
%   log c(kappa_j,beta_j) - kappa_j*r_j - beta_j*dq_j.
%
% Along the ray beta_j = e_j*kappa_j, with e_j = beta_j/kappa_j fixed at
% its unrestricted estimate, this is the function of one variable
%
%   f_j(kappa) = log c(kappa, e_j*kappa) - kappa*(r_j + e_j*dq_j),
%
% which is convex (log c is the log-partition function of an exponential
% family, hence jointly convex in (kappa,beta)) and is minimised by the
% unrestricted estimate kappa_j. Therefore, for a given m, the minimiser
% of f_j over the interval [m, restrfactor*m] is kappa_j truncated to
% that interval, and the restricted estimates are found by minimising
%
%   L(m) = sum_j n_j * f_j( min( max(kappa_j,m), restrfactor*m ) )
%
% with respect to m. L is convex in m, and its minimiser lies in the
% interval [min_j kappa_j, max_j kappa_j / restrfactor], which is the
% search interval used by fminbnd.
%
% References:
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
%   Garcia-Escudero, L.A., Gordaliza, A., Matran, C. and Mayo-Iscar, A.
%   (2008), A General Trimming Approach to Robust Cluster Analysis,
%   Annals of Statistics, Vol. 36, pp. 1324-1345.
%
%   Fritz, H., Garcia-Escudero, L.A. and Mayo-Iscar, A. (2013), A fast
%   algorithm for robust constrained clustering, Computational Statistics
%   and Data Analysis, Vol. 61, pp. 124-136.
%
% See also: tclustkent, KentMstepFS, KentlogC, restrKappaFS
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('restrKentFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    %% Two Kent groups with very different concentration.
    % The unrestricted estimates have a ratio between the kappas close to
    % 20: with restrfactor=5 the two kappas are moved towards each other
    % while the ratios beta/kappa (the shapes) do not change.
    rng(1)
    n = 300;
    Y1 = Kentrnd(n, 200, 80, [1 0 0], [0 0 1]);
    Y2 = Kentrnd(n,  10,  3, [0 1 0], [0 0 1]);
    [~,~,~,k1,b1,r1,dq1] = KentMstepFS(Y1, ones(n,1));
    [~,~,~,k2,b2,r2,dq2] = KentMstepFS(Y2, ones(n,1));
    kappa = [k1; k2]; beta = [b1; b2];
    [ko, bo] = restrKentFS(kappa, beta, [r1; r2], [dq1; dq2], [n; n], 5);
    disp(array2table([kappa ko beta bo beta./kappa bo./ko], 'VariableNames', ...
        {'kappa','kappaRestr','beta','betaRestr','shape','shapeRestr'}))
    fprintf('Ratio between restricted kappas: %.4f\n', max(ko)/min(ko));
%}
%
%{
    % restrfactor=1: common kappa.
    rng(2)
    n = 300;
    Y1 = Kentrnd(n, 100, 30, [1 0 0], [0 0 1]);
    Y2 = Kentrnd(n,  40, 15, [0 1 0], [0 0 1]);
    Y3 = Kentrnd(n,  20,  2, [0 0 1], [1 0 0]);
    [~,~,~,k1,b1,r1,dq1] = KentMstepFS(Y1, ones(n,1));
    [~,~,~,k2,b2,r2,dq2] = KentMstepFS(Y2, ones(n,1));
    [~,~,~,k3,b3,r3,dq3] = KentMstepFS(Y3, ones(n,1));
    [ko, bo] = restrKentFS([k1;k2;k3], [b1;b2;b3], [r1;r2;r3], [dq1;dq2;dq3], [n;n;n], 1);
    disp('Unrestricted kappa'); disp([k1 k2 k3])
    disp('Restricted kappa (restrfactor=1)'); disp(ko')
%}
%
%{
    % Estimates which already satisfy the constraint are not modified,
    % and an empty group (NaN) is passed through.
    kappa = [30; NaN; 15];
    beta  = [ 8; NaN;  3];
    r     = [0.96; NaN; 0.93];
    dq    = [0.02; NaN; 0.03];
    ni    = [100; 0; 100];
    [ko, bo] = restrKentFS(kappa, beta, r, dq, ni, 3);
    disp([kappa ko beta bo])
%}

%% Beginning of code

if nargin < 7 || isempty(tol)
    tol = 1e-8;
end

kappa = kappa(:); beta = beta(:); r = r(:); dq = dq(:); ni = ni(:);
kappaout = kappa;
betaout  = beta;

% Groups which take part in the restriction
active = ~isnan(kappa) & ~isnan(beta) & ni > 0 & kappa > 0;
if ~any(active)
    return
end

ka = kappa(active);
% make sure that 0 <= beta < kappa/2
ba = min(max(beta(active), 0), 0.5*ka*(1-1e-8));
ea = ba./ka;                        % shape, kept fixed
sa = r(active) + ea.*dq(active);    % sufficient statistic along the ray
na = ni(active);
betaout(active) = ba;

% Constraint already satisfied: nothing to do
if max(ka) <= restrfactor*min(ka)
    return
end

% L(m) is convex and its minimiser lies in [min(ka), max(ka)/restrfactor]
mlo = min(ka);
mhi = max(ka)/restrfactor;

    function val = Lm(m)
        kc = min(max(ka, m), restrfactor*m);
        val = 0;
        for jj = 1:numel(kc)
            val = val + na(jj)*(KentlogC(kc(jj), ea(jj)*kc(jj)) - kc(jj)*sa(jj));
        end
    end

if mhi - mlo <= tol*mlo
    mstar = mlo;
else
    mstar = fminbnd(@Lm, mlo, mhi, optimset('TolX', tol*mhi, 'Display', 'off'));
end

kc = min(max(ka, mstar), restrfactor*mstar);
kappaout(active) = kc;
betaout(active)  = ea.*kc;

end
%FScategory:CLUS-RobClaMULT
