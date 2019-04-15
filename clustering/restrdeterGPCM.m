function [lmdc]  = restrdeterGPCM(GAM, OMG, SigmaB, niini, pa)
%restrdeterGPCM applies determinat restrictions
%
%     GAM : constrained shape matrix. Matrix of size p-by-k containing in
%           column $j$, ($j=1, 2, \ldots, k$), the elements on the main
%           diagonal of shape matrix $\Gamma_j$. The elements of GAM
%           satisfy the following constraints:
%           The product of the elements of each column is equal to 1.
%           The ratio of the elements of each row is not greater than pa.shb.
%           The ratio of the elements of each column is not greater than
%           pa.shw. All the columns of matrix GAM are equal if the second
%           letter of modeltype is E. All the columns of matrix GAM are
%           equal to 1 if the second letter of modeltype is I. This matrix
%           can be constructed from routine restrshapepars
%    Omega : p-by-p-by-k 3D array containing in position j the rotation
%           matrix $\Omega_j$ for group $j$, with $j=1, 2, \ldots, k$
%   SigmaB : p-by-p-by-k array containing the k covariance matrices for the
%           k groups.
%   niini  : vector of length k containing the size of the groups.
%     pa : structure containing 3 letter character specifying modeltype,
%            number of dimensions, number of groups...
%            The fields of pars which are used in this routine are pa.p,
%            pa.k, pa.cdet
%
% 
% Output:
%
% lmdc  : row vector of length $k$ containing restricted determinants. More
%         precisely, the $j$-th element of lmdc contains $\lambda_j^{1/p}$.
%         The elements of lmdc satisfy the constraint pa.cdet in the sense that
%         $\max(lmdc)/\min(lmdc) \leq pa.cdet^{(1/p)}. In other words, the
%         ratio between the largest and the smallest determinant is not
%         greater than pa.cdet. All the elements of vector lmdc are equal
%         if modeltype is E** or if pa.cdet=1;

%% Beginning of code
% Initialize constrained determinant vector
lmd = NaN(1,pa.K);


% Inefficient code to obtain lmd
% for j=1:pa.K
%     lmd(j) = sum( diag(  diag(1./GAM(:,j)) * (OMG(:,:,j))' * SigmaB(:,:,j) * OMG(:,:,j) )) / pa.p;
% end


for j=1:pa.K
    OMGj=OMG(:,:,j);
    lmd(j) = sum( diag(OMGj' * SigmaB(:,:,j) * OMGj)./GAM(:,j) )/ pa.p;
end

% Note that ((OMG(:,:,j))' * SigmaB(:,:,j) * OMG(:,:,j) computes (\lambda_j^(1/p)*\Gamma_j) 
% where \Gamma_j is the UNCONSTRAINED shape matrix for component j. Using Luis
% Angel notation this is ( d_j^(1/p) D_j^*)
% The input of this procedure GAM, on the other hand, contains in the j-th
% column the elements on the diagonal of the CONSTRAINED shape matrix for
% component j. In Luis Angel notation the j-th column of matrix GAM contains
% d_{j1}^{***}, \ldots, d_{jp}^{***}

% lmdc = row vector containing the restricted determinants
lmdc = restreigen(lmd,niini', pa.cdet^(1/pa.p),pa.zerotol);

end
