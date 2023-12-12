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
%           the estimate of lambda. In the first iteration the estiamte of
%           \lambda_j is $|\Sigma_j|^(1/v)$, $j=1, 2, \ldots, k$ if
%           different determinants are allowed else it is a row vector of
%           ones.
%    Omega : Rotation. 3D array.
%           v-by-v-by-k 3D array containing in
%           position (:,:,j) the rotation
%           matrix $\Omega_j$ for group $j$, with $j=1, 2, \ldots, k$.
%   SigmaB : initial unconstrained covariance matrices. v-by-v-by-k array.
%            v-by-v-by-k array containing the k unconstrained covariance
%            matrices for the k groups.
%   niini  : size of the groups. Vector.
%           Row vector of length k containing the size of the groups.
%     pa : constraining parameters. Structure. Structure containing 3 letter character specifying modeltype,
%            number of dimensions, number of groups...
%            pa must contain the following fields:
%            pa.v = scalar, number of variables.
%            pa.k = scalar, number of groups.
%            pa.pars = type of Gaussian Parsimonious Clustering Model.
%               A 3 letter word in the set:
%               'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',
%               'EVI','VEI','EEI','VII','EII'
%            pa.shb = between groups shape constraint
%            pa.shw = within groups shape constraint
%            pa.zerotol = tolerance to decleare elements equal to 0.
%            pa.maxiterS = maximum number of iterations in presence of
%            varying shape matrices.
%            pa.userepmat = scalar (if =2 implicit expansion is used)
%                 Data Types - struct
%
%
%  Optional input arguments:
%
%
% Output:
%
%     GAMc : constrained shape matrix. Matrix of size v-by-k containing in
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
%   Garcia-Escudero L.A., Mayo-Iscar, A. and Riani M. (2020). Model-based
%   clustering with determinant-and-shape constraint, Statistics and
%   Computing, vol. 30, pp. 1363–1380,
%   https://link.springer.com/article/10.1007/s11222-020-09950-w
%     
%   Garcia-Escudero L.A., Mayo-Iscar, A. and Riani M. (2022). Constrained
%   parsimonious model-based clustering, Statistics and Computing, vol. 32,
%   https://doi.org/10.1007/s11222-021-10061-3
%
%
% Copyright 2008-2023.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('restrshapeGPCM')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:

%% Beginning of code
k=pa.k;
v=pa.v;
pars=pa.pars;

if strcmp(pars(2),'E') || pa.shb==1
    % In this case just restriction shw is used
    shw=pa.shw;
    
    sumnini=sum(niini);
    GAMpooled=zeros(v,v);
    
    for j=1:k
        GAMj=(niini(j)/sumnini)*  (1./lmd(j)) * (Omega(:,:,j))' * SigmaB(:,:,j) * Omega(:,:,j);
        if sum(sum(isnan(GAMj)))>0
            GAMj=zeros(v,v);
        end
        GAMpooled=GAMpooled+GAMj;
    end
    
    % Apply the constraint shw to the diagonal elements of the pooled Gamma matrix
    % GAMpooledc is a column vector of length p which contains the diagonal
    % elements of the pooled \Gamma matrix after applying constraint shw
    GAMpooledc = restreigen(diag(GAMpooled),1,shw,pa.zerotol,pa.userepmat);
    
    % Impose the constraint that the product of the elements of vector
    % GAMpooledc is equal to 1
    lmd=(prod(GAMpooledc,1)).^(1/pa.v);
    es=lmd;
    es(es==0)=1;
    GAMc=GAMpooledc./es;
    
    % replicate v-by-1 GAMpooledc vector k times
    GAMc=repmat(GAMc,1,k);
    
elseif strcmp(pars(2),'I') || pa.shw ==1
    GAMc=ones(v,k);
    
else % This is the case strcmp(pars(2),'V')
    GAM =NaN(v,k);
    shw=pa.shw;
    shb=pa.shb;
    maxiterS=pa.maxiterS;
    zerotol=pa.zerotol;
    
    for j=1:k
        Omegaj=Omega(:,:,j);
        GAM(:,j) = diag( Omegaj' * SigmaB(:,:,j) * Omegaj )/lmd(j);
    end
    
    GAMc = restrshapecore(GAM,niini,shw,shb,zerotol,maxiterS,pa.tolS,pa.sortsh,pa.userepmat);
end

end

function [GAMctrSRT]  = restrshapecore(GAMini, niini, shw, shb, zerotol, maxiterS, tolS, sortsh, userepmat)
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
% GAMini  : matrix of size v-by-k
%           GAMini contains in the first column the elements on the
%           diagonal of GAM(:,:,1). In other words the first
%           column contains the diagonal elements of matrix
%           $\Gamma_1$
%           GAMini contains in the second column the elements on the
%           diagonal of GAM(:,:,2). In other words the first
%           column contains the diagonal elements of matrix
%           $\Gamma_2$
%           .........
%   niini  : vector of length k containing the size of the groups.
%     shw  : scalar greater or equal 1. Constraint to impose inside each
%           group. For example, if shw is 3 the ratio of each column of output
%           matrix GAMc will not be greater than 3.
%     shb  : scalar greater or equal 1. Constraint to impose among groups.
%           For example, if shb is 5 the ratio of each row of (sorted) output
%           matrix GAMc will not be greater than 5.
%  zerotol : scalar. Tolerance value to declare all input values equal to 0
%           in the eigenvalues restriction routine (file restreigen.m).
% maxiterS : maximum number of iterations in the iterative procedure.
%           Scalar.
%     tolS : tolerance to use to exit the iterative procedure. Scalar. The
%           iterative procedure stops when the relative difference of
%           output matrix GAMc is smaller than itertol in two consecutive
%           iterations or when maxiterS is reached.
%
%
% Output:
%
%     GAMc : constrained shape matrix. Matrix of size v-by-k containing in
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
lamGAMc = GAMini;
% Initialize GAMc = Shape matrix constrained
[p,K] = size(GAMini);

% Apply eigenvalue restriction inside each group using constraining parameter
% shw
% The ratio of each column of matrix lamGAMc is not greater than shw
booeig=max(GAMini,[],1)./min(GAMini,[],1)>shw;
for j=1:K
    if booeig(j)==true
        lamGAMc(:,j) = restreigen(GAMini(:,j),1,shw,zerotol,userepmat);
    else
        lamGAMc(:,j) = GAMini(:,j);
    end
end


% Apply the iterative procedure to find constrained \Gamma matrix
iter=0;
% diffGAM value of the relative sum of squares of the difference between
% the element of matrix \Gamma in two consecutive iterations
diffGAM = 9999;

Ord=NaN(p,K);
GAMsor=Ord;
GAMctr=Ord;
GAMctrSRT=Ord;


while ( (diffGAM > tolS) && (iter < maxiterS) )
    iter = iter + 1;
    
    % In this stage GAM(:,j) is diag(\lambda_j^(1/p) \times \Gamma_j )
    GAM =  lamGAMc;
    GAMold=GAM(:);
    
    lmd=(prod(GAM,1)).^(1/p);
    es=lmd;
    es(es==0)=1;
    % In this stage GAM(:,j) is diag(\Gamma_j )
    % det(Gamma_j)=1
    GAM=GAM./repmat(es,p,1);
    GAM(GAM==0)=1;
    usesor=true;
    sortsh=1;
    if usesor==true
        if sortsh==1
            for j=1:K
                [GAMsor(:,j), Ord(:,j)]=sort(GAM(:,j),'ascend');
            end
        else
            GAMsor=GAM;
        end
    else
        GAMsor=GAM;
    end
    
    % Apply restriction between groups
    % The elements of each column of GAM are sorted from largest to smallest
    % The ranks of the orginal ordering of each column is store in matrix Ord
    % The ratio of each row of matrix GAMc is not greater than shb
    booeig=max(GAMsor,[],2)./min(GAMsor,[],2)>shb;
    for i=1:p
        if booeig(i)==true
            GAMctr(i,:) = restreigen(GAMsor(i,:), niini, shb, zerotol,userepmat);
        else
            GAMctr(i,:) = GAMsor(i,:);
        end
    end
    
    if usesor==true
        % Reconstruct the original order for each column
        if sortsh==1
            for j=1:K
                GAMctrSRT(Ord(:,j),j)=GAMctr(:,j);
            end
        else
            GAMctrSRT=GAMctr;
        end
    else
        GAMctrSRT=GAMctr;
    end
    
    lamGAMc=GAMctrSRT;
    
    % GAMold = old values of matrix Gamma in vectorized form
    GAMnew=lamGAMc(:);
    % diff = (new values of Gamma) - (old values of Gamma)
    diff=GAMnew-GAMold;
    % relative sum of squares of the differences
    diffGAM=diff'*diff/(GAMold'*GAMold);
    
end
end
%FScategory:CLUS-RobClaMULT
