function [Omega, Omega_]  = cpcE(lmd, SigmaB, niini, pa)

% Input p-by-p-by-K empirical covariance matrix SigmaB for the k groups
% Output 1-K vector of lambdas
% Omega_ pxp common rotation matrix
% Omega same thing as Omega_ but of size pxpxK

%% Beginning of code
p=pa.p;
K=pa.K;
sumnini=sum(niini);

Sigma_ = NaN(p,p,K);
Omega=Sigma_;
% d_=NaN(p,K);

for k=1:K
    Sigma_(:,:,k) = (1/lmd(k)) * (niini(k) /sumnini)  * SigmaB(:,:,k);
end

Sigma = sum(Sigma_,3);
[V,~]= eig(Sigma);
V=fliplr(V);
Omega_=V;
% diagSigma=diag(Sigma);
for k=1:K
    Omega(:,:,k)=V;
   % d_(:,k)=diagSigma;
end

end




