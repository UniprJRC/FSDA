function GAMc  = restrshapeGPCM(lmd, Omega, SigmaB, niini, pa)
%restrshapeGPCM produces the restricted shape matrix for the 14 GPCM
%
%
%<a href="matlab: docsearchFS('restrshapeGPCM')">Link to the help function</a>
%
%
%
% The purpose of this routine is to produce the constrained shape matrix
% $\Gamma$.
% This routine copes with the second of the 3 letters of modeltype. It
% deals with the cases in which the second letter is E, or I or V. If the
% second letter is V procedure restrshapecore is invoked and both (within
% groups) cshw, and (between groups) cshb constraints are imposed. If the
% second letter of modeltype is E just cshw is used. If the second letter
% is I, GAMc becomes a matrix of ones.
%
%
% Required input arguments:
%
%     lmd : Determinants. Vector.
%            Row vector of length k containing in the j-th position
%           $|\Sigma_j|^(1/p)$, $j=1, 2, \ldots, k$ if different
%           determinants are allowed else it is a row vector of ones.
%    Omega : Rotation. 3D array.
%           p-by-p-by-k 3D array containing in
%           position (:,:,j) the rotation
%           matrix $\Omega_j$ for group $j$, with $j=1, 2, \ldots, k$.
%   SigmaB : initial unconstrained covariance matrices. p-by-p-by-k array.
%            p-by-p-by-k array containing the k unconstrained covariance
%            matrices for the k groups.
%   niini  : size of the groups. Vector.
%           Row vector of length k containing the size of the groups.
%     pa : constraining parameters. Structure. Structure containing 3
%          letter character specifying modeltype,
%          number of dimensions, number of groups...
%          pa must contain the following fields:
%          pa.p = scalar, number of variables.
%          pa.k = scalar, number of groups.
%          pa.pars = type of Gaussian Parsimonious Clustering Model.
%               A 3 letter word in the set:
%               'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',
%               'EVI','VEI','EEI','VII','EII'
%            pa.shb = between groups shape constraint
%            pa.shw = within groups shape constraint
%            pa.zerotol = tolerance to decleare elements equal to 0.
%            pa.maxiterS = maximum number of iterations in presence of
%            varying shape matrices.
%                 Data Types - struct
%
%
%  Optional input arguments:
%
%
% Output:
%
%     GAMc : constrained shape matrix. Matrix of size p-by-k containing in
%           column j the elements on the main diagonal of shape matrix
%           $\Gamma_j$. The elements of GAMc satisfy the following
%           constraints:
%           The product of the elements of each column is equal to 1.
%           The ratio among the largest elements of each column is
%           not greater than pa.shb.
%           The ratio among the second largest elements of each column is
%           not greater than pa.shb.
%           ....
%           The ratio among the smallest elements of each column is
%           not greater than pa.shb.
%           The ratio of the elements of each column is not greater than
%           pa.shw.
%           All the columns of matrix GAMc are equal if the second
%           letter of modeltype is E. All the columns of matrix GAMc are
%           equal to 1 if the second letter of modeltype is I. This matrix
%           will be an input of routine restrdeterGPCM.m to compute
%           constrained determinants.
%
%
% See also: restrSigmaGPCM, restrdeterGPCM, restreigen, tclust
%
%
% References:
%
%   Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani M. (2019),
%   Robust parsimonious clustering models. Submitted.
%
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('restrshapeGPCM')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:

%% Beginning of code

k=pa.K;
p=pa.p;
pars=pa.pars;
if strcmp(pars(2),'E')
    % In this case just restriction shw is used
    shw=pa.shw;
    
    if sum(isnan(lmd))>0
        for j = k
            lmd(j) = (det(SigmaB(:,:,j))).^(1 /pa.p);
        end
    end
    
    sumnini=sum(niini);
    GAMpooled=zeros(p,p);
    
    for j=1:k
        GAMj=(niini(j)/sumnini)*  (1./lmd(j)) * (Omega(:,:,j))' * SigmaB(:,:,j) * Omega(:,:,j);
        if sum(sum(isnan(GAMj)))>0
            GAMj=zeros(p,p);
        end
        GAMpooled=GAMpooled+GAMj;
    end
    
    % Apply the constraint shw to the diagonal elements of the pooled Gamma matrix
    % GAMpooledc is a column vector of length p which contains the diagonal
    % elements of the pooled \Gamma matrix after applying constraint shw
    GAMpooledc = restreigen(diag(GAMpooled),1,shw,pa.zerotol);
    
    % Impose the constraint that the product of the elements of vector
    % GAMpooledc is equal to 1
    lmd=(prod(GAMpooledc,1)).^(1/pa.p);
    es=lmd;
    es(es==0)=1;
    GAMpooledc=GAMpooledc./es;
    
    % replicate p-by-1 GAMpooledc vector k times
    GAMc=repmat(GAMpooledc,1,k);
    
elseif strcmp(pars(2),'I')
    GAMc=ones(p,k);
    
else % This is the case strcmp(pars(2),'V')
    lamGAM =NaN(p,k);
    shw=pa.shw;
    shb=pa.shb;
    maxiterS=pa.maxiterS;
    zerotol=pa.zerotol;
    for j=1:k
        lamGAM(:,j) = diag( (Omega(:,:,j))' * SigmaB(:,:,j) * Omega(:,:,j) );
    end
    GAMc = restrshapecore(lamGAM,niini,shw,shb,zerotol,maxiterS,pa.tol);
end


end

function [GAMc]  = restrshapecore(lamGAM, niini, shw, shb, zerotol, maxiterS, itertol)
% restrshapecore computes constrained Gamma (shape) matrix
%
% The purpose is to find the new constrained shape matrix.
% This routine is called when modeltype(2) = V (that is in presence of
% varying shape matrices). In order to impose both shw inside each group
% and shb between groups, a number of iterations specified by input
% parameter maxiterS and a stopping condition given by itertol are necessary .
%
% Required input arguments:
%
% lamGAM  : matrix of size p-by-k
%           lamGAM contains in the first column the elements on the
%           diagonal of (lam(1)* GAM(:,:,1)). In other words the first
%           column contains the diagonal elements of matrix
%           $\lambda_1^(1/p)*\Gamma_1$
%           lamGAM contains in the second column the elements on the
%           diagonal of (lam(2)* GAM(:,:,2)). In other words the first
%           column contains the diagonal elements of matrix
%           $\lambda_2^(1/p)*\Gamma_2$
%           .........
%   niini  : vector of length k containing the size of the groups.
%     shw  : scalar greater or equal 1. Constraint to impose inside each
%           group. For example, if shw is 3 the ratio of each column of output
%           matrix GAMc will not be greater than 3.
%     shb  : scalar greater or equal 1. Constraint to impose among groups.
%           For example, if shb is 5 the ratio of each row of output
%           matrix GAMc will not be greater than 5.
%  zerotol : scalar. Tolerance value to declare all input values equal to 0
%           in the eigenvalues restriction routine (file restreigen.m).
% maxiterS : maximum number of iterations in the iterative procedure.
%           Scalar.
%  itertol : tolerance to use to exit the iterative procedure. Scalar. The
%           iterative procedure stops when the relative difference of
%           output matrix GAMc is smaller than itertol in two consecutive
%           iterations or when maxiterS is reached.
%
%
% Output:
%
%     GAMc : constrained shape matrix. Matrix of size p-by-k containing in
%           column j the elements on the main diagonal of shape matrix
%           $\Gamma_j$. The elements of GAMc satisfy the following
%           constraints:
%           The product of the elements of each column is equal to 1.
%           The ratio among the largest elements of each column is
%           not greater than pa.shb.
%           The ratio among the second largest elements of each column is
%           not greater than pa.shb.
%           ....
%           The ratio among the smallest elements of each column is
%           not greater than pa.shb.
%           The ratio of the elements of each column is not greater than
%           pa.shw.

%% Beginning of code
lamGAMc = lamGAM;
% Initialize GAMc = Shape matrix constrained
GAMc =  lamGAM;
[p,K] = size(lamGAM);
% Apply eigenvalue restriction inside each group using constraining parameter
% shw
% The ratio of each column of matrix lamGAMc is not greater than shw
for j=1:K
    lamGAMc(:,j) = restreigen(lamGAM(:,j),1,shw,zerotol);
end


% Apply the iterative procedure to find constrained \Gamma matrix
iter=0;
% diffGAM value of the relative sum of squares of the difference between
% the element of matrix \Gamma in two consecutive iterations
diffGAM = 9999;

while ( (diffGAM > itertol) && (iter < maxiterS) )
    iter = iter + 1;
    
    % In this stage GAM(:,j) is diag(\lambda_j^(1/p) \times \Gamma_j )
    GAM =  lamGAMc;
    lmd=(prod(GAM,1)).^(1/p);
    es=lmd;
    es(es==0)=1;
    % In this stage GAM(:,j) is diag(\Gamma_j )
    % det(Gamma_j)=1
    GAM=GAM./repmat(es,p,1);
    GAM(GAM==0)=1;
    
  
    % TODO TODO
    Ord=zeros(p,K);
    % GAMtmp=GAM;
    for j=1:K
        [GAM(:,j), Ord(:,j)]=sort(GAM(:,j),'desc');
    end
    
    % Apply restriction between groups
    % The elements of each column of GAM are sorted from largest to smallest
    % The ranks of the orginal ordering of each column is store in matrix ord
    % The ratio of each row of matrix GAMc is not greater than shb
    for i=1:p
        GAMc(i,:) = restreigen(GAM(i,:), niini, shb, zerotol);
    end
    
    % Reconstruct the original order for each column
    for j=1:K
        GAMc(Ord(:,j),j)=GAMc(:,j);
    end
    
    % GAMold = old values of matrix Gamma in vectorized form
    GAMold=lamGAMc(:);
    % diff = (new values of Gamma) - (old values of Gamma)
    diff=GAMc(:)-GAMold;
    % relative sum of squares of the differences
    diffGAM=diff'*diff/(GAMold'*GAMold);
    
    lamGAMc = GAMc;
end

end
%FScategory:CLUS-RobClaMULT