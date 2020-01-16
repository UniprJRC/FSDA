function [Omega, Omega2D]  = cpcV(lmdc, GAMc, Omega2D, Wk, wk, pa)
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
%           pa.shw. All the columns of matrix GAMc are equal if the second
%           letter of modeltype is E. All the columns of matrix GAMc are
%           equal to 1 if the second letter of modeltype is I.
%   Omega2D : p-by-p matrix containing the common rotation matrix.
%   SigmaB : p-by-p-by-k array containing the k covariance matrices for the
%           k groups.
%   niini  : vector of length k containing the size of the groups.
%     pa : structure containing: 3 letter character specifying modeltype,
%            number of dimensions, number of groups...
%            The fields of pars which are used in this routine are pa.p,
%            pa.k,  pa.maxiterR and pa.itertol
%
% Output:
%
%    Omega : p-by-p-k 3D array containing the updated common rotation
%               matrix replicated k times. Omega(:,:,j)=Omega2D with j=1,
%               ..., k
%
%   Omega2D : p-by-p matrix containing the updated common rotation matrix.
%
%
% References
%
% McNicholas, P.D., Browne, R.P. (2014), Estimating common principal
% components in high dimensions, Advances in Data Analysis and
% Classification, Vol. 8, pp. 217-226

% Copyright 2008-2019.
% Written by FSDA team
%$LastChangedDate::                      $: Date of the last commit


%% Beginning of code

p=pa.p;
k=pa.K;
maxiterR=pa.maxiterR;
% Wk = NaN(p,p,k);
% wk = NaN(1,k);
Fk = NaN(p, p, k);

Omega2Dini=Omega2D;

Omega = Fk;
% sumnini=sum(niini);

% Get the tolerance
itertol=pa.tol;

% Apply the iterative procedure to find Omega2D  matrix
iter=0;
% diffOMG value of the relative sum of squares of the difference between
% the element of matrix Omega2D in two consecutive iterations
diffOMG = 9999;
tmp=9999;

while ( (diffOMG > itertol) && (iter < maxiterR) )
    iter=iter+1;
    
    F=0;
    for j=1:k

        
        GAMjinv=diag(1./(GAMc(:,j)));
        F=F+(1/lmdc(j)) * GAMjinv * (Omega2D') * Wk(:,:,j)...
            - wk(j)*(lmdc(j)^(-1/p))*GAMjinv * (Omega2D');
        
    end
    
    [U,~,V] = svd(F,'econ');
    Omega2D = V*U;
    
    % Omega2Dold = old values of matrix Omega2D in vectorized form
    Omega2Dold=tmp;
    
    % Omega2new = new values of matrix Omega2D in vectorized form
    Omega2new=abs(Omega2D(:));
    % diff = (new values of Omega) - (old values of Omega)
    diff=Omega2new-Omega2Dold;
    % relative sum of squares of the differences
    % Note that (Omega2Dold'*Omega2Dold)=p
    diffOMG=diff'*diff/p;
    
    tmp = abs(Omega2new);
end
if p>2
    corm=corr(Omega2Dini,Omega2D);
    corm=abs(corm);
    [maxcor,colindmaxcor]=max(corm,[],2);
    % seqp=(p:-1:1)';
    seqp=(1:p)';
    if isequal(sort(colindmaxcor,'ascend'),seqp)
        Omega2D=Omega2D(:,colindmaxcor);
    else
        oldnewv=zeros(p,2);
        for i=1:p
            [~,indrow]=max(maxcor);
            indcol=colindmaxcor(indrow);
            oldnewv(i,:)=[indrow, indcol];
            % Put -Inf for row indrow and column indcol
            corm(indrow,:)=-Inf;
            corm(:,indcol)=-Inf;
            [maxcor,colindmaxcor]=max(corm,[],2);
        end
        oldnewv=sortrows(oldnewv,1,'ascend');
        Omega2D=Omega2D(:,oldnewv(:,2));
    end
end
%  corr(Omega2Dini,Omega2D)

% dd=1;

% Replicate Omega2D  k times
for j=1:k
    Omega(:,:,j) = Omega2D;
end
end




