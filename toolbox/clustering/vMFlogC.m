function logCp = vMFlogC(kappa,p)
%vMFlogC computes the log of the von Mises-Fisher normalizing constant Cp(kappa)
%
%<a href="matlab: docsearchFS('vMFlogC')">Link to the help page for this function</a>
%
%   vMFlogC computes, elementwise and entirely in log-space, the
%   log-normalizing constant of the von Mises-Fisher distribution on the
%   hypersphere S^(p-1):
%
%       log C_p(kappa) = (p/2-1)*log(kappa) - (p/2)*log(2*pi) - log(I_{p/2-1}(kappa))
%
%   where I_nu denotes the modified Bessel function of the first kind of
%   order nu. The logarithm of the Bessel function is computed by
%   logBesseliFS, which uses the exponentially-scaled form besseli(...,1)
%   to avoid overflow for large kappa, and an asymptotic expansion when
%   the Bessel function underflows, as it happens for large p (see the
%   section "Overflow and underflow protection" below).
%
%   Special case kappa=0 (or kappa<0, treated as 0): the vMF distribution
%   degenerates to the uniform distribution on S^(p-1) with density
%   1/omega_{p-1}, where omega_{p-1} = 2*pi^(p/2)/Gamma(p/2) is the
%   surface area of S^(p-1). The corresponding log-constant is computed
%   exactly as:
%
%       log C_p(0) = gammaln(p/2) - log(2) - (p/2)*log(pi)
%
%   This function factors out the normalizing-constant computation shared
%   by logvMFpdfFS (which needs it once per group per E-step iteration)
%   and restrKappaFS (which needs it to evaluate the vMF-consistent
%   per-group profile-likelihood loss used to restrict the concentration
%   parameters).
%
%  Required input arguments:
%
%       kappa: Concentration parameter(s). Vector or scalar. Non-negative
%              values; kappa<0 is treated as kappa=0 (uniform). Fully
%              vectorized: an entire vector of group-specific kappa values
%              can be processed in a single call.
%           p: Dimension of the ambient space. Scalar. Data lie on
%              S^(p-1), so Y has p columns. Must be an integer >=2.
%
%
%  Optional input arguments:
%
% (none beyond what is listed above; all three arguments are positional)
%
%  Output:
%
%       logCp: Log normalizing constant. Vector or scalar, same size as
%              kappa. Each entry is log C_p(kappa_i): a negative number
%              (since C_p integrates to 1 over the sphere, so C_p<1/omega
%              <1 for kappa>0, and C_p=1/omega_{p-1} for kappa=0).
%
% More About:
%
% SPECIAL CASES
%
% For p=2 (von Mises distribution on the circle), nu=0 and
%   C_2(kappa) = 1 / (2*pi*I_0(kappa)).
%
% For p=3 (spherical vMF), nu=1/2. Using I_{1/2}(k) = sqrt(2/(pi*k))*sinh(k):
%   C_3(kappa) = kappa / (4*pi*sinh(kappa)).
%
% These can be verified against logvMFpdfFS output or computed directly.
%
% OVERFLOW AND UNDERFLOW PROTECTION
%
% For large kappa (e.g. kappa>700 in double precision), I_nu(kappa)
% overflows to Inf. The scaled form besseli(nu,kappa,1)=exp(-kappa)*I_nu(kappa)
% remains in [0,1] and is then un-scaled in log-space by adding kappa back.
% For large p the opposite problem arises: the order nu=p/2-1 is large
% and I_nu(kappa) is smaller than the smallest positive double unless
% kappa is also large (for example for p=768 and kappa=50, or p=1536 and
% kappa=300). In this case besseli returns 0, its logarithm is -Inf and
% the normalizing constant would be +Inf. logBesseliFS detects this case
% and uses the uniform asymptotic expansion for large order, so that
% vMFlogC can be used with embeddings in hundreds or thousands of
% dimensions.
%
% References:
%
%   Mardia, K.V. and Jupp, P.E. (2000), Directional Statistics, Wiley,
%   Chichester. (vMF density and normalizing constant, Chapter 9)
%
%   Banerjee, A., Dhillon, I.S., Ghosh, J. and Sra, S. (2005),
%   Clustering on the Unit Hypersphere using von Mises-Fisher
%   Distributions, Journal of Machine Learning Research, Vol. 6,
%   pp. 1345-1382. (log-space evaluation of C_p(kappa))
%
% See also: logvMFpdfFS, kappaFS, restrKappaFS, logBesseliFS
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('vMFlogC')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Log normalizing constant for p=2 (von Mises) and p=3 (spherical vMF).
    % Compare against the analytic special-case formulae.
    kappa = [0, 0.5, 1, 5, 20, 100];
    logC2 = vMFlogC(kappa, 2);
    logC3 = vMFlogC(kappa, 3);
    % Analytic p=2: log C_2(k) = -log(2*pi) - log(I_0(k))
    logC2_ref = -log(2*pi) - log(besseli(0, kappa, 1)) - kappa;
    logC2_ref(kappa==0) = -log(2*pi);    % I_0(0)=1
    % Analytic p=3: log C_3(k) = log(k) - log(4*pi) - log(sinh(k))  [k>0]
    logC3_ref = nan(size(kappa));
    logC3_ref(kappa==0) = log(1/(4*pi));
    kp = kappa(kappa>0);
    logC3_ref(kappa>0) = log(kp) - log(4*pi) - log(sinh(kp));
    fprintf('p=2 max|error| = %.2e\n', max(abs(logC2 - logC2_ref)));
    fprintf('p=3 max|error| = %.2e\n', max(abs(logC3 - logC3_ref)));
%}
%
%{
    % Verify the kappa=0 formula against the sphere surface area.
    % For kappa=0, C_p(0) = 1/omega_{p-1} = Gamma(p/2)/(2*pi^(p/2)).
    for p = [2, 3, 4, 5, 10, 20]
        logC0 = vMFlogC(0, p);
        omega  = 2*pi^(p/2) / gamma(p/2);    % surface area of S^(p-1)
        logC0_ref = -log(omega);
        fprintf('p=%2d: vMFlogC=%.6f, -log(omega)=%.6f, err=%.1e\n', ...
            p, logC0, logC0_ref, abs(logC0 - logC0_ref));
    end
%}
%
%{
    % Overflow protection: large kappa does not produce Inf or NaN.
    % Direct evaluation of I_{p/2-1}(kappa) overflows for large kappa,
    % but the scaled form in vMFlogC remains well-behaved.
    p = 3;
    kappa = [1e2, 1e3, 1e4, 1e5];
    logCp = vMFlogC(kappa, p);
    fprintf('kappa=%7.0f: logCp = %.4f\n', [kappa; logCp]);
    % Approximate for large kappa: log C_p(k) ~ (p/2-1)*log(k) - (p/2)*log(2*pi)
    %   - k + 0.5*log(2*pi/k)   [from the large-k Bessel asymptotics]
    logCp_approx = (p/2-1)*log(kappa) - (p/2)*log(2*pi) ...
                   - kappa + 0.5*log(2*pi./kappa);
    fprintf('  large-kappa approx: %.4f\n', logCp_approx);
%}
%
%{
    %% High dimension: the constant is finite when the Bessel function underflows.
    % p=1536 (dimension of common document embeddings). The value is
    % checked by numerical integration: the density of t=mu'*y is
    % C_p(kappa)*omega_{p-2}*exp(kappa*t)*(1-t^2)^((p-3)/2), where
    % omega_{p-2} is the surface area of S^(p-2), and it integrates to 1.
    p = 1536; kappa = 300;
    logC = vMFlogC(kappa, p);
    logomega = log(2) + ((p-1)/2)*log(pi) - gammaln((p-1)/2);
    f = @(t) exp(logC + logomega + kappa*t + ((p-3)/2)*log1p(-t.^2));
    fprintf('log C_p(kappa) = %.4f\n', logC);
    fprintf('Integral of the density of mu''*y: %.8f\n', integral(f, 0, 0.5));
%}

%% Beginning of code

nu = p/2 - 1;
logCp = zeros(size(kappa));

pos = kappa>0;
if any(pos)
    kp = kappa(pos);
    % log(I_nu(kappa)): reliable also when I_nu(kappa) underflows (large
    % p) or overflows (large kappa)
    logIbessel = logBesseliFS(nu, kp);
    logCp(pos) = (p/2 - 1).*log(kp) - (p/2)*log(2*pi) - logIbessel;
end
if any(~pos)
    logCp(~pos) = -( log(2) + (p/2)*log(pi) - gammaln(p/2) );
end

end
%FScategory:CLUS-RobClaMULT
