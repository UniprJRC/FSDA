function [lmd]  = restrdeterpars(GAM, OMG, SigmaB, niini, pa, zero_tol)

% Find best determinant vector
% 

%% Beginning of code

if nargin<5
    zero_tol=1e-16;
end

% autovalues = lmd;
autovalues = NaN(1,pa.K);
for k=1:pa.K
    autovalues(k) = sum(diag(  diag(1./GAM(:,k)) * (OMG(:,:,k))' * SigmaB(:,:,k) * OMG(:,:,k) )) / pa.p;
end
lmd = restreigen(autovalues,niini', pa.cdet^(1/pa.p),zero_tol);

end
