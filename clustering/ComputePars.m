function [li,di,const1]=ComputePars(v, k, Pi, Mu, S, S05, Sinv, detS)
%ComputePars computes parameters needed for computing overlap
%
%  Required input arguments:
%
%       v  : scalar, dimensionality (number of variables)
%       k  : scalar, number of components (groups)
%       Pi : vector of length k containing mixing proportions
%       Mu : matrix of size k-by-v containing centroids of the k groups
%       S  : 3D array of size v-by-v-by-k containing covariance matrices
%           S= Sigma_1, ... , \Sigma_k
%
%  Optional input arguments:
%
%      S05 : 3D array of size v-by-v-by-k containing matrices 
%           \Sigma_1^{0.5} ... , \Sigma_k^{0.5}
%
%     Sinv : 3D array of size v-by-v-by-k containing matrices 
%           \Sigma_1^{-1} ... , \Sigma_k^{-1}
%
%     detS : vector of length k containing determinants of the covariance
%            matrices of the k groups
%            detS= (|\Sigma_1|, |\Sigma_2|, ..., |\Sigma_k|)   
%    
%           Remark: if optional input arguments are supplied then the
%           fifth input argument S is unnecessary and is ignored
%
%  Output:
%
%       li : 3D array of size k-by-k-by-v. li(i,j,:) is the vector which
%            contains the eigenvalues of matrix Sji=\Sigma_j|i
%            where \Sigma_j|i  = \Sigma_i^0.5 \Sigma_j^{-1} \Sigma_i^0.5.
%            Remark: li(i,j,:)-1 are the coefficients of the linear
%            combinations of non central chi^2 distribution which is used
%            to compute the probability of overlapping
%            The non centrality parameter  of the non central chi^2 distribution
%            associated with li(i,j,l)-1 is (li(i,j,l) di(i,j,l)/(li(i,j,l)-1))^2
%
%       di : 3D array of size k-by-k-by-p
%            di(i,j,l)= \gamma_l' \Sigma_i^{-0.5}(\mu_i-\mu_j)
%            l=1, 2, ...,v, where \gamma_l is the l-th eigenvector coming
%            from the spectral decomposition of matrix \Sigma_j|i  . Vector
%            di(i,j,:) contains the coefficients of N(0,1) (called W_l in
%            equation 2.2 of Maitra and Melnikov (2010), JCGS)
%
%    const1: k x k matrix the  element i,j with i=1,..., k and j=1, ..., k
%            is equal to log((Pi(j)/Pi(i))^2 * detS(i)/detS(j))
%            where detS(i) and detS(j) are, respectively, the
%            determinants of the covariance matrices of groups i and j. See
%            equation (2.2) of Maitra and Melnikov (2010), JCGS
%            Note that const1(j,i)=-const1(i,j). The elements on the
%            diagonal of matrix const1 are set to 0 (in other words they
%            are not computed).
%
%            The misclassification probability \omega_{j|i}=
%               + \sum_{l=1}^p (li(i,j,l)-1) U_l                         (1)
%               + 2 \sum_{l=1}^p di(i,j,l) \W_l                          (2)
%                            \leq
%               + \sum_{l=1}^p \frac{li(i,j,l) *di(i,j,l)^2}{li(i,j,l-1) (3)
%               - \sum_{l=1}^p di(i,j,l)^2                               (4)
%                + log const1(i,j)
%              where U_l are non central \chi^2 r.v. with 1 degree of freedom and
%              non centrality parameter (li(i,j,l)^2 di(i,j,l)^2/(li(i,j,l)-1)^2
%              To be precise, sums in (1) and (3) are for l:li(i,j,l) is different from 1
%              Similarly, sums in (2) and (4) are for l:li(i,j,l) is = 1
%
% Copyright 2008-2019.
%
%$LastChangedDate::                      $: Date of the last commit

% Examples:

%{
    rng(123,'twister');
    k=4; % Number of groups
    p=5;    % Number of dimensions
    Pi=[0.1 0.2 0.4 0.3]; % mixing proportions
    % Mu matrix of means of size k-by-p; (each row is a distinct centroid)
    Mu=randn(k,p);
    % Groups 2 and 3  is far from the other groups
    Mu(2:3,:)=Mu(2:3,:)+10;
    %Mu(3,:)=Mu(3,:)+2;
    %Mu(4,:)=Mu(4,:)+3;
    % S= 3D array of dimension p-by-p-by-k containing covariance matrices of
    % the groups
    S=zeros(p,p,k);
    for j=1:k
        S(:,:,j)=eye(p);
    end

    [li, di, const1]=ComputePars(p, k, Pi, Mu, S);

    asympt = 0;
    c = 1;
    fixcl=zeros(k,1);
    tol=1e-8;
    lim=1e+07;
    [OmegaMap,  BarOmega, MaxOmega, rcMax] = ...
        GetOmegaMap(c, p, k, li, di, const1, fixcl, tol, lim, asympt);
    disp('Omegamap= k-by-k matrix which contains misclassification probabilities')
    disp(OmegaMap);
    disp('Average overlap')
    disp(BarOmega)
    disp('Maximum overlap')
    disp(MaxOmega)
    disp('Groups with maximum overlap')
    disp(rcMax)

%}

%% Beginning of code

if nargin==5
    % S05 = p-by-p-by-k which will contain in the j-th third dimension
    % \Sigma_j^{0.5} j=1, 2, ..., k
    S05=zeros(v,v,k);
    % Sinv = p-by-p-by-k which will contain in the j-th third dimension
    % \Sigma_j^{-1} j=1, 2, ..., k
    Sinv=S05;
    % detS = vector of length k which will contain the determinant of the
    % covariance matrices
    detS=zeros(k,1);
    
    for kk=1:k
        
        % Sk = covariance matrix of group k
        Sk=S(:,:,kk);
        % The line below is added to ensure numerical stability in
        % computing the eigenvalues and eigenvectors, in the case Si has
        % entries very close to 0.
        Sk(abs(Sk)<1e-10)=0;
        
        % ViSk= matrix containing the eigevectors of Sk (covariance matrix of
        % group k)
        % eigSk = diagonal matrix containing the eigenvelues of Sk (covariance
        % matrix of group k)
        [ViSk,eigSk]=eig(Sk);
        eigSk=real(eigSk);
        
        % detS(kk) = determinant of covariance matrix of group k
        detS(kk) = prod(diag(eigSk));
        
        % SingvalSk = matrix containing singular values (sqrt of the eigenvalues) of
        % covariance matrix of group k
        SingvalSk=eigSk^0.5;
        
        % S05(:,:,kk) contains matrix \Sigma_k^{0.5}
        S05(:,:,kk)=ViSk*SingvalSk*(ViSk');
        
        % eigSinvk = eigenvalues of matrix \Sigma_k^{-1}
        eigSinvk=diag(diag(eigSk).^-1);
        
        % Sinv(:,:,kk) contains matrix \Sigma_k^{-1}
        Sinv(:,:,kk)=ViSk*eigSinvk*(ViSk');
    end
end

li=zeros(k,k,v);
di=li;
const1=zeros(k,k);

for ii=1:k-1
    for jj=ii+1:k
        % S05i = \Sigma^{0.5}_i
        S05i=S05(:,:,ii);
        % Sinvj =  \Sigma^{-1}_j
        Sinvj=Sinv(:,:,jj);
        % Sinvi = \Sigma^{-1}_i
        Sinvi=Sinv(:,:,ii);
        
        m1=(Mu(ii,:)-Mu(jj,:))';
        m2=-m1;
        
        % Sji is matrix \Sigma^{0.5}_i \Sigma^{-1}_j \Sigma^{0.5}_i
        % This matrix is called \Sigma_{j|i}
        Sji=S05i*Sinvj*(S05i);
        
        % The line below is added to ensure numerical stability in
        % computing the eigenvalues and eigenvectors, in the case Si has
        % entries very close to 0.
        Sji(abs(Sji)<1e-10)=0;
        
        % Spectral decomposition of matrix Sji = \Sigma_{j|i}
        % ViSji= \Gamma_j|i (matrix of eigenvectors of Sji)
        % EigSji = \Lambda_j|i (matrix of eigenvalues of Sji)
        [ViSji,EigSji]=eig(Sji);
        
        % Store in third dimension of array li
        % the p eigenvalues of matrix Sji =\Sigma_{j|i}
        li(ii,jj,:)=diag(EigSji);
        
        % Sinvi = \Sigma_i^{-1}
        % S05i =\Sigma_i^{0.5}
        % S05invi= \Sigma_i^{-0.5} = \Sigma_i^{-1} * \Sigma_i^{0.5}
        S05invi=Sinvi*S05i;
        
        % ViSji = (delta_1, ..., delta_l, ... \delta_k)
        %  di(ii,jj,l) = \gamma_l \Sigma_i^{-0.5} (\mu_i -\mu_j) l=1, 2, ...,p
        % See Theorem 1 of Maitra and Melnykov (2010), JCGS
        di(ii,jj,:)= (ViSji')*(S05invi*m1);
        
        % Lines below find di(jj,ii,:) and li(jj,ii,:)
        
        % S05j = \Sigma^{0.5}_j
        S05j=S05(:,:,jj);
        % Sij = \Sigma^{0.5}_j   \Sigma^{-1}_i \Sigma^{0.5}_j
        % This matrix is called \Sigma_{i|j}
        Sij=S05j*Sinvi*S05j;
        
        % The line below is added to ensure numerical stability in
        % computing the eigenvalues and eigenvectors, in the case Si has
        % entries very close to 0.
        Sij(abs(Sij)<1e-10)=0;
        % Compute eigenvectors and eigenvalues of matrix Sij=\Sigma_{i|j}
        [ViSij,EigSij]=eig(Sij);
        
        % Store in third dimension of array li
        % the p eigenvalues of matrix Sij =\Sigma_{i|j}
        li(jj,ii,:) = diag(EigSij);
        
        % Sinvj = \Sigma_j^{-1}
        % S05j =\Sigma_j^{0.5}
        % S05invj= \Sigma_j^{-0.5} = \Sigma_j^{-1} * \Sigma_j^{0.5}
        S05invj=Sinvj*S05j;
        di(jj,ii,:)= (ViSij')*(S05invj*m2);
        
        const1(ii,jj) = log((Pi(jj)/Pi(ii))^2 * detS(ii)/detS(jj));
        const1(jj,ii) = -const1(ii,jj);
        
    end
end

end