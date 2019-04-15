function [Omega, Omega2D]  = cpcV(lmdc, GAMc, Omega2D, SigmaB, niini, pa)
%cpcV computes updated common rotation matrix when shapes are different
%
%  This routine is called when the parametrization is *VE that is when
%  variable shape is assumed but equal rotation is imposed. This routine is
%  based on the algorithm described in McNicholas and Browne (2014)
%
% Required input arguments:
%
%
%     lmdc  : row vector of length $k$ containing restricted determinants. More
%             precisely, the $j$-th element of lmdc contains $\lambda_j^{1/p}$.
%             The elements of lmdc satisfy the constraint pa.cdet in the sense that
%             $\max(lmdc)/\min(lmdc) \leq pa.cdet^{(1/p)}. In other words, the
%             ratio between the largest and the smallest determinant is not
%             greater than pa.cdet. All the elements of vector lmdc are equal
%             if modeltype is E** or if pa.cdet=1;
%     GAMc : constrained shape matrix. Matrix of size p-by-k containing in
%           column j the elements on the main diagonal of shape matrix
%           $\Gamma_j$. The elements of GAMc satisfy the following
%           constraints:
%           The product of the elements of each column is equal to 1.
%           The ratio of the elements of each row is not greater than pa.shb.
%           The ratio of the elements of each column is not greater than
%           pa.shw. All the columns of matrix GAM are equal if the second
%           letter of modeltype is E. All the columns of matrix GAM are
%           equal to 1 if the second letter of modeltype is I. 
%   Omega2D : p-by-p matrix containing the common rotation matrix.
%   SigmaB : p-by-p-by-k array containing the k covariance matrices for the
%           k groups.
%   niini  : vector of length k containing the size of the groups.
%     pa : structure containing: 3 letter character specifying modeltype,
%            number of dimensions, number of groups...
%            The fields of pars which are used in this routine are pa.p,
%            pa.k and pa.maxiterR
%
% Output:
%
%    Omega : p-by-p-k 3D array containing the updated common rotation
%               matrix replicated k times. Omega(:,:,j)=Omega2D with j=1,
%               ..., k
%   Omega2D : p-by-p matrix containing the updated common rotation matrix.
%
%
% References 
%
% McNicholas, P.D., Browne, R.P. (2014), Estimating common principal
% components in high dimensions, Advances in Data Analysis and
% Classification, Vol. 8, pp. 217-226

%% Beginning of code
p=pa.p;
k=pa.K;
maxiterR=pa.maxiterR;
Wk = NaN(p,p,k);
wk = NaN(1,k);
Fk = NaN(p, p, k);
Omega = Fk;
sumnini=sum(niini);
for i=1:maxiterR
    for j=1:k
        Wk(:,:,j)  = (niini(j) /sumnini)  * SigmaB(:,:,j);
        wk(j) = max(eig(Wk(:,:,j)));
        Fk(:,:,j) = (1/lmdc(j)) * diag(1./(GAMc(:,j))) * (Omega2D') * Wk(:,:,j)...
            - wk(j)*(lmdc(j)^(-1/ p))*diag(1./GAMc(:,j)) * (Omega2D');
    end
    F = sum(Fk,3);
    [U,~,V] = svd(F);
    Omega2D = V*U;
end
for j=1:k
    Omega(:,:,j) = Omega2D;
end
end




