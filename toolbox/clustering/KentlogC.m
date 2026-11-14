function logc = KentlogC(kappa, beta, mterms)
%KentlogC computes the log of the Kent distribution normalizing constant log c(kappa,beta)
%
%<a href="matlab: docsearchFS('KentlogC')">Link to the help page for this function</a>
%
%   KentlogC evaluates, in log-space, the normalizing constant of the Kent
%   (Fisher-Bingham FB5) distribution on S^2 (the unit sphere in R^3):
%
%       c(kappa,beta) = 2*pi * sum_{j=0}^{inf}
%                           Gamma(j+1/2)/Gamma(j+1) * beta^(2j)
%                           * (kappa/2)^(-2j-1/2) * I_{2j+1/2}(kappa)
%
%   as given in Kent (1982, eq. 3) and Mardia et al. (2021, eq. 5), where
%   I_nu is the modified Bessel function of the first kind. Each term is
%   evaluated in log-space using the exponentially scaled Bessel function
%   besseli(nu,kappa,1), so that no overflow occurs for large kappa, and
%   the terms are added with the log-sum-exp device.
%
%   Special cases:
%     beta=0:  the Kent density reduces to the vMF density with parameter
%              kappa, and
%              log c(kappa,0) = log(4*pi) + log(sinh(kappa)) - log(kappa)
%              = -vMFlogC(kappa,3).
%     kappa=0: (with beta=0) uniform distribution on S^2, c = 4*pi.
%
%   NOTE: the Kent distribution is defined ONLY on S^2 (p=3). This
%   function does not generalise to higher-dimensional spheres (unlike the
%   vMF functions in this toolbox).
%
%  Required input arguments:
%
%       kappa: Concentration parameter. Scalar >= 0.
%        beta: Ovalness parameter. Scalar >= 0. The Kent density is
%              unimodal for 0 <= beta < kappa/2. beta=0 gives the
%              isotropic (vMF) special case.
%
%  Optional input arguments:
%
%      mterms: Number of terms of the series. Scalar positive integer or
%              empty. If it is empty or missing (default) the number of
%              terms is chosen automatically: terms are added until they
%              become negligible with respect to the partial sum. If
%              mterms is given, exactly mterms terms (j=0,...,mterms-1)
%              are used; note that a fixed small number of terms (e.g.
%              the value 20 suggested by Mardia et al., 2021) is not
%              enough when kappa is large and beta is close to kappa/2.
%                Example - 30
%                Data Types - double
%
%  Output:
%
%        logc: log c(kappa,beta). Scalar. The Kent log-density is
%              -logc + kappa*(gamma3'*x) + beta*[(gamma1'*x)^2-(gamma2'*x)^2].
%
% More About:
%
% CONVERGENCE OF THE SERIES
%
% The ratio between two consecutive terms is
%
%   t_j/t_{j-1} = (2j-1)/(2j) * (2*beta/kappa)^2 * I_{2j+1/2}(kappa)/I_{2j-3/2}(kappa)
%
% (the order of the Bessel function increases by 2 from one term to the
% next). All factors are smaller than 1 when beta<kappa/2, so the terms
% decrease monotonically, but the decrease is slow when 2*beta/kappa is
% close to 1 and kappa is large: for this reason the default is an
% adaptive number of terms. The automatic rule stops when the last term
% is smaller than eps*(1-(2*beta/kappa)^2) times the partial sum (the
% second factor bounds the geometric tail), or after 5000 terms.
%
% VERY LARGE KAPPA
%
% For kappa>1e8 the Bessel functions cannot be computed accurately and
% the high-concentration approximation of Kent (1982)
%   c(kappa,beta) ~ 2*pi*exp(kappa)/sqrt((kappa-2*beta)*(kappa+2*beta))
% is returned.
%
% References:
%
%   Kent, J.T. (1982), The Fisher-Bingham distribution on the sphere,
%   J. Roy. Statist. Soc., B, Vol. 44, pp. 71-80.
%
%   Mardia, K.V., Barber, S., Burdett, P.M., Kent, J.T. and Hamelryck, T.
%   (2021), Mixture models for spherical data with applications to protein
%   bioinformatics, arXiv:2104.13140.
%
% See also: logKentpdfFS, KentMstepFS, restrKentFS, tclustkent, Kentrnd
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('KentlogC')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Compare KentlogC with the analytic vMF special case (beta=0).
    % For beta=0, log c(kappa,0) = log(4*pi) + log(sinh(kappa)) - log(kappa).
    kappa = [1 5 20 100];
    beta  = zeros(size(kappa));
    logc_kent = arrayfun(@(k,b) KentlogC(k,b), kappa, beta);
    logc_vMF  = log(4*pi) + log(sinh(kappa)) - log(kappa);
    fprintf('kappa=%6.1f: Kent series=%.6f, vMF analytic=%.6f, diff=%.2e\n', ...
        [kappa; logc_kent; logc_vMF; abs(logc_kent-logc_vMF)]);
%}
%
%{
    %% Compare KentlogC with numerical integration of the Kent kernel.
    % In the canonical frame, with t=cos(colatitude) and longitude phi,
    % c(kappa,beta) = int_0^{2pi} int_{-1}^{1}
    %                 exp(kappa*t + beta*(1-t^2)*cos(2*phi)) dt dphi.
    % The factor exp(-kappa) is taken out to avoid overflow.
    kb = [5 2; 20 9; 30 10; 100 45];
    for i = 1:size(kb,1)
        kappa = kb(i,1); beta = kb(i,2);
        f = @(t,phi) exp(kappa*(t-1) + beta*(1-t.^2).*cos(2*phi));
        logc_num = log(integral2(f,-1,1,0,2*pi,'AbsTol',1e-13,'RelTol',1e-11)) + kappa;
        fprintf('kappa=%5.1f beta=%5.1f: series=%.6f, integral=%.6f\n', ...
            kappa, beta, KentlogC(kappa,beta), logc_num);
    end
%}
%
%{
    % Case kappa=0: the constant is 4 pi and logc equals log(4 pi).
    logc0 = KentlogC(0, 0);
    fprintf('KentlogC(0,0) = %.6f,  log(4*pi) = %.6f\n', logc0, log(4*pi));
%}
%
%{
    % Fixed versus automatic number of terms.
    % With kappa large and beta close to kappa/2 a small fixed number of
    % terms is not enough.
    kappa = 1000; beta = 480;
    for m = [10 20 40 80 160]
        fprintf('mterms=%3d: logc = %.10f\n', m, KentlogC(kappa, beta, m));
    end
    fprintf('automatic : logc = %.10f\n', KentlogC(kappa, beta));
%}
%
%{
    % logc increases with beta (kappa fixed).
    kappa = 30;
    betas = [0 2 5 10 14];   % all < kappa/2 = 15
    for b = betas
        fprintf('beta=%4.1f: logc = %.4f\n', b, KentlogC(kappa,b));
    end
%}

%% Beginning of code

if nargin < 3
    mterms = [];
end

% Degenerate case: kappa=0 -> uniform on S^2, c=4*pi
if kappa <= 0
    logc = log(4*pi);
    return
end

% Extremely large kappa: besseli loses all its significant digits, and
% the leading term of the asymptotic expansion of Kent (1982) is used
%   c(kappa,beta) ~ 2*pi*exp(kappa)/sqrt((kappa-2*beta)*(kappa+2*beta))
if kappa > 1e8
    logc = log(2*pi) + kappa ...
        - 0.5*log(max(kappa-2*beta, eps*kappa)) - 0.5*log(kappa+2*beta);
    return
end

adaptive = isempty(mterms);
if beta <= 0
    % only the term j=0 is different from 0
    adaptive = false;
    chunk = 1;
elseif adaptive
    chunk = 50;
else
    chunk = mterms;
end
maxterms = 5000;

if beta > 0
    log2bk = log(2*beta/kappa);
    % stopping rule of the adaptive version (see More About)
    logtol = log(eps) + log(max(1 - (2*beta/kappa)^2, 1e-3));
end

% log of term j:
%   gammaln(j+1/2) - gammaln(j+1) + 2j*log(beta) - (2j+1/2)*log(kappa/2)
%   + log I_{2j+1/2}(kappa)
% = gammaln(j+1/2) - gammaln(j+1) + 2j*log(2*beta/kappa) - 0.5*log(kappa/2)
%   + log(besseli(2j+1/2,kappa,1)) + kappa
logS = -Inf;
j0 = 0;
while true
    j = (j0:j0+chunk-1)';
    lt = gammaln(j+0.5) - gammaln(j+1) - 0.5*log(kappa/2) ...
        + log(besseli(2*j+0.5, kappa, 1)) + kappa;
    if beta > 0
        lt = lt + 2*j*log2bk;
    end
    % log-sum-exp of the partial sum and of the new terms
    mx = max(max(lt), logS);
    logS = mx + log(sum(exp(lt-mx)) + exp(logS-mx));
    j0 = j0 + chunk;
    if ~adaptive || j0 >= maxterms || lt(end) < logS + logtol
        break
    end
end

logc = log(2*pi) + logS;

end
%FScategory:CLUS-RobClaMULT
