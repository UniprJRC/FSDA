function [lmdc]  = restrdeterGPCM(GAM, OMG, SigmaB, niini, pa)
%restrdeterGPCM applies determinat restrictions for the 14 GPCM
%
%
%<a href="matlab: docsearchFS('restrdeterGPCM')">Link to the help function</a>
%
%
%  This routine applies the constraints to the determinants using the
%  specification contained in field pa.cdet of input structure pa.
%
%
%  Required input arguments:
%
%     GAM : constrained shape matrix. 2D array.
%           Matrix of size p-by-k containing in
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
%           Data Types - double
%    OMG  : costrained rotation array. 3D array. p-by-p-by-k 3D array
%           containing in position (:,:,j) the rotation
%           matrix $\Omega_j$ for group $j$, with $j=1, 2, \ldots, k$
%           Data Types - double
%   SigmaB : initial unconstrained covariance matrices. p-by-p-by-k array.
%            p-by-p-by-k array containing the k unconstrained covariance
%            matrices for the k groups.
%   niini  : size of the groups. Vector.  
%           Row vector of length k containing the size of the groups.
%           Data Types - double
%     pa : constraining parameters. Structure. Structure containing 3 letter character specifying modeltype,
%            number of dimensions, number of groups...
%            pa must contain the following fields: 
%            pa.v = scalar, number of variables.
%            pa.k = scalar, number of groups.
%            pa.cdet = determinants constraint
%           Data Types - double
% 
%
%  Optional input arguments:
%
%
% Output:
%
% lmdc  : restricted determinants. Vector. 
%         Row vector of length $k$ containing restricted determinants. More
%         precisely, the $j$-th element of lmdc contains $\lambda_j^{1/p}$.
%         The elements of lmdc satisfy the constraint pa.cdet in the sense
%         that $\max(lmdc) / \min(lmdc) \leq pa.cdet^{(1/p)}$. In other words,
%         the ratio between the largest and the smallest determinant is not
%         greater than pa.cdet. All the elements of vector lmdc are equal
%         if modeltype is E** or if pa.cdet=1;
%
%
%
% See also: restrSigmaGPCM, restrdeterGPCM, restreigen, tclust
%
%
% References:
%
%   Garcia-Escudero L.A., Mayo-Iscar, A. and Riani M. (2020). Model-based
%   clustering with determinant-and-shape constraint, Statistics and
%   Computing, vol. 30, pp. 1363–1380,
%   https://link.springer.com/article/10.1007/s11222-020-09950-w
%     
%   Garcia-Escudero L.A., Mayo-Iscar, A. and Riani M. (2022). Constrained
%   parsimonious model-based clustering, Statistics and Computing, vol. 32,
%   https://doi.org/10.1007/s11222-021-10061-3
% 
% Copyright 2008-2025.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('restrdeterGPCM')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:

%% Beginning of code
% Initialize constrained determinant vector
lmd = NaN(1,pa.k);


% Inefficient code to obtain lmd
% for j=1:pa.k
%     lmd(j) = sum( diag(  diag(1./GAM(:,j)) * (OMG(:,:,j))' * SigmaB(:,:,j) * OMG(:,:,j) )) / pa.v;
% end

% Fancy way which avoids the loop below (it works from MATLAB 2020b)
% MM=pagemtimes(pagemtimes(OMG,'ctranspose',SigmaB,'none'),OMG);
% Lr=reshape(MM,pa.v*pa.v,[]);
% D=reshape(Lr(1:(pa.v+1):end,:),pa.v,pa.k);
% lmd=sum(D./GAM,1)/pa.v;
for j=1:pa.k
    OMGj=OMG(:,:,j);
    lmd(j) = sum( diag(OMGj' * SigmaB(:,:,j) * OMGj)./GAM(:,j) )/ pa.v;
end


% Note that ((OMG(:,:,j))' * SigmaB(:,:,j) * OMG(:,:,j) computes
% (\lambda_j^(1/p)*\Gamma_j) where \Gamma_j is the UNCONSTRAINED shape
% matrix for component j based on updated Omega matrix. Using Luis Angel
% notation this is ( d_j^(1/p) D_j^*)
% GAM, the input of this procedure, on the other hand, contains in the j-th
% column the elements on the diagonal of the CONSTRAINED shape matrix for
% component j. In Luis Angel notation the j-th column of matrix GAM contains
% d_{j1}^{***}, \ldots, d_{jp}^{***}

% lmdc = row vector containing the restricted determinants
% Make sure niini is a column vector
niini=niini(:);
if max(lmd)/min(lmd)>pa.cdet^(1/pa.v)
    lmdc = restreigen(lmd,niini, pa.cdet^(1/pa.v),pa.zerotol,pa.userepmat);
else
    lmdc=lmd;
end

end

%FScategory:CLUS-RobClaMULT