function [Sigma]  = restrSigmaGPCM(SigmaB, niini, pa)
%restrpars computes updated common rotation matrix when shapes are equal
%
%
%<a href="matlab: docsearchFS('restrSigmaGPCM')">Link to the help function</a>
%
%
%  This routine applies the constraints to the covariance matrices using the
%  specifications contained in input structure pa.
%
%
% Required input arguments:
%
%   SigmaB : initial unconstrained covariance matrices. p-by-p-by-k array.
%            p-by-p-by-k array containing the k covariance matrices for the
%            k groups.
% 
%   niini  : sizes of the groups. Vector. Row vector of length k containing
%           the size of the groups.
%
%      pa  : Constraints to apply and model specification. Structure.
%            Structure containing the following fields:
%             pa.pars= type of Gaussian Parsimonious Clustering Model. Character.
%               A 3 letter word in the set:
%               'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',
%               'EVI','VEI','EEI','VII','EII'
%             pa.cdet = scalar in the interval [1 Inf) which specifies the 
%               the restriction which has to be applied to the determinants.
%               If pa.cdet=1 all determinants are forced to be equal. 
%               See section More About for additional details.
%             pa.shw = scalar in the interval [1 Inf) which specifies the 
%               the restriction which has to be applied to the elements of
%               the shape matrices inside each group. If pa.shw=1 all diagonal
%               elements of the shape matrix of cluster j (with j=1, ...,
%               k) will be equal.
%             pa.shb = scalar in the interval [1 Inf) which specifies the 
%               the restriction which has to be applied to the elements of
%               the shape matrices across each group. 
%             pa.maxiterS = positive integer which specifies the maximum
%               number of iterations to obtain the restricted shape matrix.
%               This parameter is used by routine restrshapepars. The
%               default value of pa.maxiterS is 5.
%             pa.maxiterR = positive integer which specifies the maximum
%               number of iterations to obtain the common rotation matrix
%               in presence of varying shape.
%               This parameter is used by routine cpcV. The
%               default value of pa.maxiterR is 20.
%          pa.maxiterDSR = positive integer which specifies the maximum
%               number of iterations to obtain the requested restricted
%               determinants, shape matrices and rotation. For all
%               parametrizations  pa.maxiterDSR is set to 1 apart from for
%               the specifications 'VVE', 'EVE' and 'VEE'. The default
%               value of pa.maxiterDSR is 20.
%           pa.tol=tolerance to use to exit the iterative procedures. Scalar. The
%               iterative procedures stops when the relative difference of
%               a certain output matrix is smaller than itertol in two consecutive
%               iterations. The default value of pa.tol is 1e-12.
%      pa.zerotol = tolerance value to declare all input values equal to 0
%               in the eigenvalues restriction routine (file restreigen.m)
%               or in the final reconstruction of covariance matrices. 
%               The default value of zerotol is 1e-10.
%           pa.k  = the number of groups.
%           pa.p  = the number of variables.
%
%  Optional input arguments:
%
%
% Output:
% 
%
%             Sigma  : constrained covariance matrices. p-by-p-by-k array.
%                     p-by-p-by-k array containing the k covariance
%                     matrices for the k groups.
%
% More About:
% The notation for the eigen-decomposition of the
% component covariance matrices is as follows
% 
% \[
% \Sigma_j= \lambda_j^{1/p} \Omega_j \Gamma_j \Omega_j'  \qquad j=1, 2, \ldots, k
% \]
% The dimension of matrices $\Omega_j$ and $\Lambda_j$ is $p\times p$.
% 
% $c_{det}=$ scalar, constraint associated with the determinants.
% 
% $c_{shw}=$ scalar, constraint inside each group of the shape matrix.
% 
% $c_{shb}=$ scalar, constraint among groups of the shape matrix.
% 
% Note that if you impose equal volumes $c_{det}=1$. Similarly, if you impose a spherical shape $c_{shw}=1$.
%
% We also denote with
% 
%   [1] $\Sigma$ the 3D array of size $p\times p \times k$ containing the
%   empirical covariance matrices of the $k$ groups, before applying the
%   constraints coming from the 14 parametrizations. In the code $\Sigma$
%   is called $SigmaB$. The $j$-th slice of this 3D array of size $p\times
%   p$ is denoted with symbol $\hat \Sigma_j$.
%   [2] $\Omega$ the 3D array of size $p\times p \times k$ containing the
%   rotation matrices of the $k$ groups. In the code $\Omega$  is called
%   $OMG$. The $j$-th slice of this 3D array of size $p\times p$ is called
%   $\hat \Omega_j$.
%   [3] $\Gamma$ the $p\times k$ matrix containing in column $j$, with
%   $j=1, 2, \ldots, k$, the diagonal elements of matrix $\Gamma_j$ (shape
%   matrix of group j). In the code matrix $\Gamma$ is called GAM.
%   [4] niini the vector of length $k$ containing the number of units
%   (weights) associated to each group. 
%   [5]  $\lambda$ = the vector of length $p$ containing in the $j$-th
%   position $\lambda_j^{1/p}=|\Sigma_j|^{1/p}$. In the code  vector
%   $\lambda$ is called $lmd$.
%	
% 	
%
% See also restrshape, restrdeter, restreigen, tclust
%
%
% References:
%
%   Garcia-Escudero, L.A., Mayo-Iscar, A. and Riani M. (2019),
%   A General Trimming Approach to Robust Cluster Analysis. Submitted.
%
%
% Copyright 2008-2018.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('restrSigmaGPCM')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit

% Examples:

%
%{


%}


%% Beginning of code

% SigmaB = v-times-v-times-k = empirical covariance matrix
% pa = structure containing modeltype, number of iterations .....
Sigma=SigmaB;
K=length(niini);
lmd=NaN(1,K);

% OMG = initialize 3D array containing rotation matrices
OMG=zeros(size(Sigma));

% Tolerance associated to the maximum of the elements which have to be
% constrained. If the maximum  is smaller than zerotol restreigen
% procedures returns in output what has been given in input. For example,
% if the elements which have to be constrained are the eigenvalues of the
% covariance matrices and the max of the eigenvalues is smaller than
% zerotol it means that all n points are concentrated in k points and there
% is a perfect fit therefore no further changes on the eigenvalues is
% required.
zero_tol=pa.zerotol;

% pa.pars = character vector with three letters specifying the type of the
% 14 contraints (i.e. EEE, CVVV, EVE, ...)
pars=pa.pars;


% VEV VEE EVE VVE require iterations because in these cases common rotation
% matrix is updated in each step.
% All the other specification do not
if strcmp(pars,'EEE') || strcmp(pars,'VVV') || strcmp(pars,'EVV') ...
        || strcmp(pars,'EEV') || strcmp(pars(3),'I')
    pa.maxiterDSR = 1;
end

% If equal determinants are imposed pa.cdet=1
if strcmp(pars(1),'E')
    pa.cdet = 1;
end

% If OMG is identity, shape restriction parameter within groups is set to 1
if strcmp(pars(2),'I')
    pa.shw = 1;
end

%% Initialization part
if strcmp(pars(3),'E')
    % In the common principal components case it is necessary to find
    % initial values for OMG, GAM and lmd
    
    % Find initial values of lmd (unconstrained determinants) and OMG (rotation)
    [lmd, OMG]  = initR(SigmaB, niini, pa);
    % Find initial constrained shape matrix GAM
    % pa.shw and pa.shb constraining parameters are used
    [GAM] =restrshapeGPCM(lmd, OMG, SigmaB, niini, pa);
    % Find initial constrained determinants (lmd vector)
    lmd =restrdeterGPCM(GAM, OMG, SigmaB, niini, pa);
    
elseif  strcmp(pars(3),'V')
    % Find initial (and final value for OMG)
    for j=1:K
        [V,~]= eig(SigmaB(:,:,j));
        V=fliplr(V);
        % D=fliplr(D);
        OMG(:,:,j)=V;
    end
    
else % The remaning case is when **I
    % Find initial (and final value for OMG).
    % in this case OMG is a 3D arry contaning identity matrices
    eyep=eye(pa.p);
    for j=1:K
        OMG(:,:,j)=eyep;
    end
end

% End of initialization

%% Beginning of iterative process
for i=1:pa.maxiterDSR
    
    % In the **E case (except for the case EEE) it is necessary to update in each step of the
    % iterative procedure OMG 
    
    if strcmp(pars,'VVE') || strcmp(pars,'EVE')
        % Variable shape: update OMG (rotation)
        % parameter pa.maxiterR is used here
        [OMG]  = cpcV(lmd, GAM, OMG(:,:,1), SigmaB, niini, pa);
        
    elseif strcmp(pars,'VEE') || strcmp(pars,'EEE') % TODO check if EEE is necessary here because given that all lmd are equal the order of the eigenvectors can change
        % Equal shape: update OMG (rotation)
        % old eiggiroe
        [OMG]  = cpcE(lmd, SigmaB, niini, pa);
    else
        % In all the other cases OMG is not updated
    end
    
    % Find new value of GAM (shape matrix) using updated OMG
    [GAM] =restrshapeGPCM(lmd, OMG, SigmaB, niini, pa);
    % Find new value of lmd (determinants) starting from constrained shape
    % GAM and updated OMG
    [lmd]=restrdeterGPCM(GAM, OMG, SigmaB, niini, pa);
end

% Check if all is well
codenonZero = max(max(GAM)) > zero_tol;

if  codenonZero
    % Reconstruct the cov matrices using final values of lmd, OMG and GAM
    for j=1:K
        Sigma(:,:,j) = lmd(j) * OMG(:,:,j)* diag(GAM(:,j))* (OMG(:,:,j)');
    end
else
    % In this case final Sigma is the identity matrix replicated k times
    eyep=eye(p);
    for j=1:K
        Sigma(:,:,j) = eyep;
    end
end


end