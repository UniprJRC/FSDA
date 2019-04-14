function [Omega, Omega_]  = cpcV(lmd, LMD, Omega_, SigmaB, niini, pa)

% Input p-by-p-by-K empirical covariance matrix SigmaB for the k groups
% Output 1-K vector of lambdas
% Omega_ pxp common rotation matrix
% Omega same thing as Omega_ but of size pxpxK

%% Beginning of code
p=pa.p;
K=pa.K;
maxiterR=pa.maxiterR;
Wk = NaN(p,p,K);
wk = NaN(1,K);
Fk = NaN(p, p, K);
Omega = Fk;
sumnini=sum(niini);
for i=1:maxiterR
    for k=1:K
        Wk(:,:,k)  = (niini(k) /sumnini)  * SigmaB(:,:,k);
        wk(k) = max(eig(Wk(:,:,k)));
        Fk(:,:,k) = (1/lmd(k)) * diag(1./(LMD(:,k))) * (Omega_') * Wk(:,:,k)...
            - wk(k)*(lmd(k)^(-1/ p))*diag(1./LMD(:,k)) * (Omega_');
    end
    F = sum(Fk,3);
    [U,~,V] = svd(F);
    Omega_ = V*U;
end
   for k=1:K
        Omega(:,:,k) = Omega_;
    end
end




