function [Omega, Omega2D]  = cpcE(lmdc, SigmaB, niini, pa)
%cpcE computes updated common rotation matrix when shapes are equal
%
%  This routine is called when the parametrization is *EE that is when
%  equal shape and equal rotation are imposed. 
%
%
% Required input arguments:
%
%     lmdc  : row vector of length $k$ containing restricted determinants. More
%             precisely, the $j$-th element of lmdc contains $\lambda_j^{1/p}$.
%             The elements of lmdc satisfy the constraint pa.cdet in the sense that
%             $\max(lmdc)/\min(lmdc) \leq pa.cdet^{(1/p)}. In other words, the
%             ratio between the largest and the smallest determinant is not
%             greater than pa.cdet. All the elements of vector lmdc are equal
%             if modeltype is E** or if pa.cdet=1;
%   SigmaB : p-by-p-by-k array containing the k covariance matrices for the
%           k groups.
%   niini  : vector of length k containing the size of the groups.
%     pa : structure containing: 3 letter character specifying modeltype,
%            number of dimensions, number of groups...
%            The fields of pars which are used in this routine are pa.p,
%            and pa.k a
%
%
% Output:
%
%    Omega : p-by-p-k 3D array containing the updated common rotation
%               matrix replicated k times. Omega(:,:,j)=Omega2D with j=1,
%               ..., k
%   Omega2D : p-by-p matrix containing the updated common rotation matrix.
%
% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit

%% Beginning of code

p=pa.p;
k=pa.K;
sumnini=sum(niini);
% Inefficient way of obtaining Sigma
% Sigma_ = NaN(p,p,k);
% Omega=Sigma_;
% for j=1:k
%     Sigma_(:,:,j) = (1/lmdc(j)) * (niini(j) /sumnini)  * SigmaB(:,:,j);
% end
% Sigma = sum(Sigma_,3);

Omega=NaN(p,p,k);
% Sigma is OMG*GAM*OMG' pooled
Sigma=zeros(p,p);
for j=1:k
    Sigma = Sigma + (1/lmdc(j)) * (niini(j) /sumnini)  * SigmaB(:,:,j);
end

[V,~]= eig(Sigma);
V=fliplr(V);
Omega2D=V;
% diagSigma=diag(Sigma);

% Create updated Omega
for j=1:k
    Omega(:,:,j)=V;
end

end




