function logI = logBesseliFS(nu, x)
%logBesseliFS computes the logarithm of the modified Bessel function of the first kind
%
%<a href="matlab: docsearchFS('logBesseliFS')">Link to the help page for this function</a>
%
%   logBesseliFS returns log(I_nu(x)), where I_nu is the modified Bessel
%   function of the first kind of order nu>=0, for x>=0. The value is
%   reliable also when I_nu(x) cannot be represented in double precision,
%   that is when it overflows (large x) or underflows (large nu). The
%   second case is the relevant one for directional data in high
%   dimension: the normalizing constant of the von Mises-Fisher
%   distribution on S^(p-1) contains I_{p/2-1}(kappa), and for p in the
%   hundreds or thousands (e.g. document embeddings) and moderate kappa
%   this quantity is smaller than the smallest positive double, so that
%   besseli returns 0 and its logarithm is -Inf.
%
%   logBesseliFS is used by vMFlogC (normalizing constant) and by kappaFS
%   (ratio of Bessel functions A_p(kappa)).
%
%  Required input arguments:
%
%           nu: Order. Non-negative scalar, or array of the same size as
%               x.
%            x: Argument. Array of non-negative values, or scalar.
%
%  Optional input arguments:
%
%  Output:
%
%         logI: log(I_nu(x)). Array of the same size as x (or as nu if x
%               is a scalar). It is -Inf when x=0 and nu>0.
%
% More About:
%
% The value is first computed as log(besseli(nu,x,1))+x, where
% besseli(nu,x,1)=exp(-x)*I_nu(x) is the exponentially scaled function,
% which does not overflow. When this is not finite (besseli underflows to
% 0, or it cannot be computed because the argument is too large) the
% following asymptotic expansions are used.
%
% For nu>=1, the uniform asymptotic expansion for large order (Debye;
% Abramowitz and Stegun, 1964, 9.7.7): with z=x/nu, s=sqrt(1+z^2), t=1/s
% and eta=s+log(z/(1+s)),
%
%   log I_nu(x) ~ nu*eta - log(2*pi*nu)/2 - log(s)/2
%                 + log(1 + u1(t)/nu + u2(t)/nu^2 + u3(t)/nu^3),
%
% where u1, u2, u3 are polynomials in t. The error of the expansion
% decreases with nu, and it is below 1e-10 in the cases in which besseli
% underflows (which require nu larger than about 100).
%
% For nu<1 (the expansion is needed only for very large x), the expansion
% for large argument (Abramowitz and Stegun, 1964, 9.7.1):
%
%   log I_nu(x) ~ x - log(2*pi*x)/2 + log(1 - (4*nu^2-1)/(8*x)).
%
% References:
%
%   Abramowitz, M. and Stegun, I.A. (1964), Handbook of Mathematical
%   Functions, National Bureau of Standards, Washington.
%
% See also: vMFlogC, kappaFS, logvMFpdfFS
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('logBesseliFS')">Link to the help page for this function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%
% Examples:
%
%{
    % Agreement with besseli when besseli can be computed.
    nu = 2.5;
    x = [0.1 1 10 100 500];
    disp([log(besseli(nu,x)); logBesseliFS(nu,x)])
%}
%
%{
    %% Large order: besseli underflows, logBesseliFS does not.
    % Order nu=p/2-1 for embeddings of dimension p=1536, with kappa=300.
    p = 1536; kappa = 300;
    nu = p/2 - 1;
    fprintf('log(besseli)  = %g\n', log(besseli(nu,kappa,1)) + kappa);
    fprintf('logBesseliFS  = %.6f\n', logBesseliFS(nu,kappa));
%}
%
%{
    % Large argument: besseli overflows without the scaling.
    nu = 0.5; x = 1e4;
    % exact value for nu=1/2: I(x)=sqrt(2/(pi*x))*sinh(x)
    exact = 0.5*log(2/(pi*x)) + x - log(2);
    fprintf('logBesseliFS = %.8f, exact = %.8f\n', logBesseliFS(nu,x), exact);
%}

%% Beginning of code

% nu and x of the same size
if isscalar(nu) && ~isscalar(x)
    nu = repmat(nu, size(x));
elseif isscalar(x) && ~isscalar(nu)
    x = repmat(x, size(nu));
end

logI = log(besseli(nu, x, 1)) + x;

% Values which besseli cannot represent
bad = ~isfinite(logI) & x > 0 & ~isnan(x) & ~isnan(nu);
if any(bad(:))
    nub = nu(bad);
    xb  = x(bad);
    lb  = zeros(size(xb));

    % large order: uniform (Debye) expansion
    big = nub >= 1;
    if any(big)
        n = nub(big);
        z = xb(big)./n;
        s = sqrt(1 + z.^2);
        t = 1./s;
        eta = s + log(z./(1 + s));
        u1 = (3*t - 5*t.^3)/24;
        u2 = (81*t.^2 - 462*t.^4 + 385*t.^6)/1152;
        u3 = (30375*t.^3 - 369603*t.^5 + 765765*t.^7 - 425425*t.^9)/414720;
        lb(big) = n.*eta - 0.5*log(2*pi*n) - 0.5*log(s) ...
            + log1p(u1./n + u2./n.^2 + u3./n.^3);
    end

    % small order and very large argument
    if any(~big)
        n = nub(~big);
        xs = xb(~big);
        lb(~big) = xs - 0.5*log(2*pi*xs) + log1p(-(4*n.^2 - 1)./(8*xs));
    end

    logI(bad) = lb;
end

end
%FScategory:UTISTAT
