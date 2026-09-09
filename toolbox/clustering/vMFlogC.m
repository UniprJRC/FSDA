function logCp = vMFlogC(kappa,p)
%vMFlogC computes the log of the von Mises-Fisher normalizing constant Cp(kappa)
%
%<a href="matlab: docsearchFS('vMFlogC')">Link to the help function</a>
%
%   vMFlogC computes, elementwise,
%
%       log Cp(kappa) = (p/2-1)*log(kappa) - (p/2)*log(2*pi) - log(I_(p/2-1)(kappa))
%
%   in log-space using the exponentially-scaled Bessel function, exactly
%   as done inside logvMFpdfFS (this function factors that computation
%   out so it can also be used by restrKappaFS, which needs the same
%   normalizing constant to build the profile-likelihood loss used to
%   restrict the concentration parameters).
%
%  Required input arguments:
%
%       kappa: Concentration parameter(s). Vector or scalar. kappa>=0.
%           p: Dimension of the ambient space. Scalar. Data lie on S^(p-1).
%
%  Output:
%
%       logCp: log Cp(kappa). Vector or scalar, same size as kappa.
%
% See also: logvMFpdfFS, kappaFS, restrKappaFS
%
% Copyright 2008-2025.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('vMFlogC')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

nu = p/2 - 1;
logCp = zeros(size(kappa));

pos = kappa>0;
if any(pos)
    kp = kappa(pos);
    logIbessel = log(besseli(nu, kp, 1)) + kp;
    logCp(pos) = (p/2 - 1).*log(kp) - (p/2)*log(2*pi) - logIbessel;
end
if any(~pos)
    logCp(~pos) = -( log(2) + (p/2)*log(pi) - gammaln(p/2) );
end

end
%FScategory:CLUS-RobClaMULT
