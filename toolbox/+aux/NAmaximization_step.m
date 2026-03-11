function [mu_new, Sigma_new] = NAmaximization_step(T1, T2, w)
% NAmaximization_step  M-step: update mu and Sigma from expected stats.
%
% Required input arguments:
%
%   T1 - p x 1 weighted sum of expected x
%   T2 - p x p weighted sum of expected xx'
%   w  - n x 1 weight vector (or a scalar sum-of-weights)
%
%  Optional input arguments:
%
%  Output:
%
% mu_new = T1 / sum(w)
% Sigma_new = T2 / sum(w) - mu_new * mu_new'
%
% Copyright 2008-2026.
% Written by FSDA team
%

%% Beginning of code
if numel(w) > 1
    s = sum(w);
else
    s = w;
end
mu_new = T1 / s;
Sigma_new = (T2 / s - (mu_new * mu_new'))*(s/(s-1));

% to ensure symmetry numerical issues -> symmetrize
Sigma_new = (Sigma_new + Sigma_new') / 2;

end