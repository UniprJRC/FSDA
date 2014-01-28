%computes parameters needed for computing overlap
%  INPUT PARAMETERS
%  * p  - dimensionality
%  * k  - number of components
%  * Pi - mixing proportions
%  * Mu - mean vectors
%  * S  - covariance matrices
%  OUTPUT PARAMETERS
%  * li, di, const1 - parameters needed for computing overlap (see theory of method)
%  */

function [li,di,const1]=ComputePars(p,k,Pi,Mu,S)

Sh=zeros(p,p,k);
detS=zeros(k,1);
Sinv=Sh;
li=zeros(k,k,p);
di=li;
const1=zeros(k,k);

for kk=1:k
    
    Ga=S(:,:,kk);
    % The line below is added to ensure numerical stability in
    % computing the eigenvalues and eigenvectors, in the case Si has
    % entries very close to 0.
    Ga(abs(Ga)<1e-10)=0;
    
    [VGa,eigGa]=eig(Ga);
    
    detS(kk) = prod(diag(eigGa));
    
    L=eigGa^0.5;
    
    % Sh contains matrices \Sigma_i^{0.5}
    Sh(:,:,kk)=VGa*L*(VGa');
    
    % Sinv contains matrices \Sigma_i^{-1}
    % L = eigenvalues of matrix \sigma_i^{-1}
    L=diag(diag(eigGa).^-1);
    
    Sinv(:,:,kk)=VGa*L*(VGa');
end


for ii=1:k-1
    for jj=ii+1:k
        % Ga = \Sigma^{0.5}_i
        Ga=Sh(:,:,ii);
        % L =  \Sigma^{-1}_j
        L=Sinv(:,:,jj);
        % L2 = \Sigma^{-1}_i
        L2=Sinv(:,:,ii);
        
        %for kk=1:p
        %    m1(kk) = Mu(ii,kk) - Mu(jj,kk);
        %    m2(kk) = -m1(kk);
        %end
        
        m1=(Mu(ii,:)-Mu(jj,:))';
        m2=-m1;
        
        % Si is matrix \Sigma^{0.5}_i \Sigma^{0.5}_j \Sigma^{0.5}_i
        Si=Ga*L*(Ga);
        
        % the line below is added to ensure numerical stability in
        % computing the eigenvalues and eigenvectors, in the case Si has
        % entries very close to 0.
        Si(abs(Si)<1e-10)=0;
        [Vi,Eigi]=eig(Si);
        
        % li(i,j) are the eigenvalues of
        % of matrix Si
        
        li(ii,jj,:)=diag(Eigi);
        
        % L2 = \Sigma_i^{-1}
        % Ga =\Sigma_i^{0.5}
        % Ga2= \Sigma_i^{-0.5}
        Ga2=L2*Ga;
        % Eig = \Sigma_i^{-0.5} *(\mu_i-\mu_j)
        Eig=Ga2*m1;
        
        % Ga2= Transpose of the matrix which contains the eigenvectors
        % of matrix Si
        Ga2=Vi';
        m1=Ga2*Eig;
        
        %multiply(L2, p, p, Ga, p, p, Ga2);
        % matxvec(Ga2, p, p, m1, p, Eig);
        % tA(Si, p, p, Ga2);
        %matxvec(Ga2, p, p, Eig, p, m1);
        
        % store the \delta_l which are an ingredient to find the non
        % centrality parameters
        di(ii,jj,:)= m1;
        
        % Remark: final non centrality parameters will be given by
        %  ........
        
        % Ga2 = \Sigma^{0.5}_j
        Ga2=Sh(:,:,jj);
        % Si = \Sigma^{0.5}_j   \Sigma^{-1}_i \Sigma^{0.5}_j
        Si=Ga2*L2*Ga2;
        % cpy1(Sh, j, p, p, Ga2);
        % XAXt(Ga2, p, L2, Si);
        
        % The line below is added to ensure numerical stability in
        % computing the eigenvalues and eigenvectors, in the case Si has
        % entries very close to 0.
        Si(abs(Si)<1e-10)=0;
        [Vi,Eigi]=eig(Si);
        
        %	EigValDec(p, Eig, Si, &dtmt);
        li(jj,ii,:) = diag(Eigi);
        
        Ga=L*Ga2;
        Eig=Ga*m2;
        Ga=Vi';
        m2= Ga*Eig;
        %multiply(L, p, p, Ga2, p, p, Ga);
        %matxvec(Ga, p, p, m2, p, Eig);
        % tA(Si, p, p, Ga);
        % matxvec(Ga, p, p, Eig, p, m2);
        
        di(jj,ii,:)= m2;
        
        const1(ii,jj) = log((Pi(jj)/Pi(ii))^2 * detS(ii)/detS(jj));
        const1(jj,ii) = -const1(ii,jj);
        
    end
end

end
