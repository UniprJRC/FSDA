function [lmd, Omega, Omega2D]  = initR(SigmaB, niini, pa)
%initR finds initial common rotation matrix
%
% This procedures is called when **E that is when a common rotation matrix
% is imposed. The main purpose of this function is to find the an initial
% estimate of the common rotation matrix.
%
%
% Required input arguments:
%
%
%   SigmaB : p-by-p-by-k array containing the k covariance matrices for the
%           k groups.
%   niini  : vector of length k containing the size of the groups.
%     pars : structure containing 3 letter character specifying modeltype,
%            number of dimensions, number of groups...
%            The fields of pars which are used in this routine are pa.p,
%            pa.k and pa.pars
%
% Output:
%
%     lmd : row vector of length k containing in the j-th position
%           $|\Sigma_j|^(1/p)$, $j=1, 2, \ldots, k$ if different
%           determinants are allowed else it is a row vector of ones.
%   Omega2D : p-by-p matrix containing the eigenvectors of pooled matrix
%            $\sum_{j=1}^k \frac{n_j}{n} \frac{1}{|\Sigma_j|^(1/p)}
%            \Sigma_j$. The first column is associated with the largest
%            eigenvalue .... This is the initial common rotation matrix.
%    Omega : p-by-p-k 3D array containing in position j Omega2D.
%            This is the common rotation matrix replicated k times.

% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit


%% Beginning of code

p=pa.p;
k=pa.K;
pars=pa.pars;

% Initialize lmd
lmd = NaN(1,k);

if strcmp(pars(1),'V')
    for j=1:k
        lmd(j) = (det(SigmaB(:,:,j))) ^ (1 / p);
    end
else
    lmd = ones(1,k);
end

sumniini=sum(niini);
% Sigma = NaN(p, p, K);
% for j=1:K
%     Sigma(:,:,j) = (niini(j) / sumniini) * (1 /lmd(j) ) * SigmaB(:,:,j);
% end
% Sw = pooled estimate of \Omega \Gamma \Omega^T (rescaled pooled
% within-group covariance matrix)
% Sw = sum(Sigma,3);
Sw=zeros(p,p);
for j=1:k
    Sw = Sw + (niini(j) / sumniini) * (1 /lmd(j) ) * SigmaB(:,:,j);
end

maxsigma_=max(max(isnan(Sw)));

% In case of missing or infinite values common rotation matrix is forced to
% be the identity matrix
if isnan(maxsigma_) || isinf(maxsigma_)
    Omega2D = eye(p);
else
    % The common rotation matrix is formed by the eigenvectors of the
    % pooled within-group covariance matrix Sw
    [V,eigunsorted]= eig(Sw);
    % Sort eigenvalues from largest to smallest and reorder the columns
    % of the matrix of eigenvectors accordingly
    [~,ordeig]=sort(diag(eigunsorted),'desc');
    V=V(:,ordeig);
    Omega2D=V;
end

% Omega = array of size p-by-p-by-k
% containing k replicates of matrix V
Omega=repmat(Omega2D,1,1,k);
end


