function [OmegaMap, BarOmega, MaxOmega, rcMax]=GetOmegaMap(c, v, k, li, di, const1, fix, tol, lim, asympt)
%GetOmegaMap calculates the map of misclassificaton between groups
%
%  Required input arguments:
%
%       c  : scalar, inflation parameter for covariance matrices. In the
%            case of heterogeneous clusters scalar c is used to correct
%            li(i,j) and const1(i,j). More precisely when hom=0 and 
%            fix(i)=0 di(i,j,:)=di(i,j,:)/c^0.5  
%                   (because di(i,j,:)=\Gamma *(c*\Sigma_i)^{-0.5}(\mu_i-\mu_j)
%                    where Gamma is the matrix of eigenvectors of
%                    \Sigma_j|i)
%                   if fix(j)=1  li(i,j,:)=c li(i,j,:), const1(i,j)=const1(i,j)+v log(c) 
%                   (because the eigenvalues of matrix 
%                   \Sigma_j|i = (c^0.5 \Sigma_i) \Sigma_j^-1 (c^0.5 \Sigma_i) are multiplied by c
%                   Similarly log |c\Sigma_i|= c*log(v) + log |\Sigma_i|
%                   if fix(j)=0  li and const1 are not changed because
%                   \Sigma_j|i remains unaltered
%       v  : dimensionality (number of variables)
%       k  : number of components (groups)
%       li : 3D array of size k-by-k-by-p
%       di : 3D array of size k-by-k-by-p
%    const1: k x k matrix
%           REMARK: li, di and const1 are the parameters needed for computing
%           overlap 
%      fix : vector of length k containing zeros or ones
%          if fix(j) =1 cluster j does not participate to inflation or
%          deflation. If fix=zeros(k,1) all clusters participate in inflation/deflation
%          This parameter is used just if heterogeneous clusters are used
%      tol : scalar. Error bound for overlap computation default is 1e-06
%      lim : maximum number of integration terms default is 1e06.
%               REMARK: Optional parameters tol and lim will be used by
%               function ncx2mixtcdf.m which computes the cdf of a linear
%               combination of non central chi2 r.v.. This is the
%               probability of overlapping
%   asympt : flag for regular or asymptotic overlap. If asypmt ==1 formula
%            of asymptotic overlap is used (see p. 359, paragraph starting
%            with step 3 of Maitra and Melnikov, 2010, JCGS). In this case
%            the misclassification probability \omega_j|i (that is with
%            respect to the centroid and covariance matrix of group i) is
%            given by 
%            \sum_{l=1}^p (\lambda_l-1) U_l \leq log \pi_j^2 |\Sigma_j| / \pi_i^2 |\Sigma_i| 
%            where U_l are central chi2 with one degree of freedom.
%            Note that parameter asympt is set to 1 in the preliminary
%            phase just to check whether the average or the maximum overlap
%            can be reached.
%
%  Output:
%
%    OmegaMap : k-by-k matrix containing map of misclassification
%               probabilities. More precisely, OmegaMap(i,j) is the
%               probability that group i overlaps with group j
%               (i ~= j)=1, 2, ..., k
%    BarOmega : scalar associated with average overlap.
%               BarOmega is computed as (sum(sum(OmegaMap))-k)/(0.5*k(k-1))
%    MaxOmega : scalar associated with maximum overlap. MaxOmega is the
%               maximum of OmegaMap(i,j)+OmegaMap(j,i)
%               (i ~= j)=1, 2, ..., k
%       rcMax : column vector of length equal to 2 containing the indexes
%               associated with the pair of components producing the
%               highest overlap (largest off diagonal element of matrix
%               OmegaMap)
%
% Copyright 2008-2014.
% Written by FSDA team
%
%{
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
    disp('Omegamap= k-by-k matrix which will contain misclassification probabilities')
    disp(OmegaMap);
    disp('Average overlap')
    disp(BarOmega)
    disp('Maximum overlap')
    disp(MaxOmega)
    disp('Groups with maximum overlap')
    disp(rcMax)
%}

%% Beginning of code

% Omegamap= k-by-k matrix which will contain misclassification probabilities
OmegaMap=eye(k);
rcMax=zeros(2,1);

TotalOmega = 0.0;
MaxOmega = 0.0;

df=ones(v,1);


ii = 1;
jj = 2;


% check if clusters are homogeneous 
hom = 1;
for kk=1:v
    
    if li(1,2,kk) ~= li(2,1,kk)
        hom = 0;
        break
    end
end


if hom == 1   % homogeneous clusters 

    if asympt == 0
        
        while ii < k
            
            Di =  squeeze(di(ii,jj,:)) / sqrt(c);
            Cnst1 = sum(Di.^2);
            
            % t is the value in which the cdf of the mixture of
            % non central chi2 must be evaluated. 
            t = const1(ii,jj) - Cnst1;
            sigma = 2 * sqrt(Cnst1);
            
            OmegaMap(ii,jj)=normcdf(t,0,sigma);
            % This is equivalent to 
            %{
                % coef = coefficients of the linear combination of non
                % central chi2 distributions (in the case of homogeneous
                % clusters coef is a vector of 0s)
                % df = degrees of freedom of mixture of non central chi2
                coef = zeros(p,1);
                ncp = zeros(p,1);
                OmegaMap(ii,jj)=ncx2mixtcdf(t, df, coef,ncp,'sigma',sigma,'lim',lim,'tol',tol);
            %}
            % see beginning of section 2.1 
            
            % t   :         scalar, value at which the cdf must be evaluated
            % df  :         vector of length k containing the degrees of freedom of the
            %               k non central chi2 distributions
            % coef  :         vector of length k containing the coefficients of the linear combinations
            %               of the k non central chi2 distributions
            %  ncp :       vector of length k containing the k non centrality parameters
            %               of the k non central chi2 distributions
            
            Di =  squeeze(di(jj,ii,:)) / sqrt(c);
            Cnst1 = sum(Di.^2);
            
            
            
            t = const1(jj,ii) - Cnst1;
            sigma = 2 * sqrt(Cnst1);
            
             OmegaMap(jj,ii)=normcdf(t,0,sigma);
            % This is equivalent to 
            %{
                % coef = coefficients of the linear combination of non
                % central chi2 distributions (in the case of homogeneous
                % clusters coef is a vector of 0s)
                % df = degrees of freedom of mixture of non central chi2
                coef = zeros(p,1);
                ncp = zeros(p,1);
                OmegaMap(jj,ii) = ncx2mixtcdf(t, df, coef,ncp,'sigma',sigma,'lim',lim,'tol',tol);
            %}
            
            OmegaOverlap = OmegaMap(ii,jj) + OmegaMap(jj,ii);
            TotalOmega = TotalOmega + OmegaOverlap;
            
            if OmegaOverlap > MaxOmega
                MaxOmega = OmegaOverlap;
                rcMax(1) = ii;
                rcMax(2) = jj;
            end
            
            
            if (jj < k) 
                jj = jj + 1;
            else
                ii = ii + 1;
                jj = ii + 1;
            end
            
        end
        
    end
    
    
    
    if asympt == 1
        
        while ii < k
            
            if const1(ii,jj) > 0
                OmegaMap(ii,jj) = 1;
                OmegaMap(jj,ii) = 0;
            end
            
            if const1(ii,jj) < 0
                OmegaMap(ii,jj) = 0;
                OmegaMap(jj,ii) = 1;
            end
            
            if const1(ii,jj) == 0
                OmegaMap(ii,jj) = 0.5;
                OmegaMap(jj,ii) = 0.5;
            end
            
            
            OmegaOverlap = OmegaMap(ii,jj) + OmegaMap(jj,ii);
            TotalOmega = TotalOmega + OmegaOverlap;
            
            if OmegaOverlap > MaxOmega
                MaxOmega = OmegaOverlap;
                rcMax(1) = ii;
                rcMax(2) = jj;
            end
            if (jj < k) %(k - 1)
                jj = jj + 1;
            else
                ii = ii + 1;
                jj = ii + 1;
            end
            
        end
        
    end
end


% hom =0 ==> non homogeneous (non equal) covariance matrices
% that is heterogeneous clusters 
if hom == 0
    
    sigma = 0.0;
    
    if (asympt == 0)
        
        while ii < k
            
            if fix(ii) == 1
                Di = squeeze(di(ii,jj,:));
                
                if fix(jj) == 1
                    Li = squeeze(li(ii,jj,:));
                    Cnst1 = const1(ii,jj);
                else
                    Li = squeeze(li(ii,jj,:))/ c;
                    Cnst1 = const1(ii,jj) - v * log(c);
                end
                
            else
                
                
                Di = squeeze(di(ii,jj,:)) / sqrt(c);
                
                if fix(jj) == 1
                    
                    Li = c * squeeze(li(ii,jj,:));
                    
                    Cnst1 = const1(ii,jj) + v * log(c);
                else
                    Li = squeeze(li(ii,jj,:));
                    
                    Cnst1 = const1(ii,jj);
                end
                
            end
            
            
            
            coef = Li - 1;
            ldprod = Li.*Di;
            const2 = (ldprod.*Di)./ coef;
            % ncp = non centrality parameters
            % The l-th element of ncl i=1, 2, ..., p is equal to
            % \lambda_l^2 \delta_l^2 /(\lambda_l-1)^2 
            ncp = (ldprod./coef).^2;
            % sum(const2)= \sum_{l=1}^p \lambda_l \delta_l^2 /(\lambda_l-1) 
            t = sum(const2)+ Cnst1;
            
            OmegaMap(ii,jj)=ncx2mixtcdf(t, df, coef, ncp,'sigma',sigma,'lim',lim,'tol',tol);

            if fix(jj) == 1
                
                
                Di = squeeze(di(jj,ii,:));
                
                if fix(ii) == 1
                    Li = squeeze(li(jj,ii,:));
                    Cnst1 = const1(jj,ii);
                else
                    Li = squeeze(li(jj,ii,:)) / c;
                    Cnst1 = const1(jj,ii) - v * log(c);
                end
                
            else
                
                Di = squeeze(di(jj,ii,:)) / sqrt(c);
                
                
                if fix(ii) == 1
                    Li = c * squeeze(li(jj,ii,:));
                    
                    Cnst1 = const1(jj,ii) + v * log(c);
                else
                    Li = squeeze(li(jj,ii,:));
                    Cnst1 = const1(jj,ii);
                end
                
            end
            
            coef = Li - 1.0;
            ldprod = Li.*Di;
            const2 = (ldprod.*Di)./ coef;
            ncp = (ldprod./coef).^2;
            t = sum(const2)+ Cnst1;
            
            OmegaMap(jj,ii)=ncx2mixtcdf(t, df, coef,ncp,'sigma',sigma,'lim',lim,'tol',tol);
            
            OmegaOverlap = OmegaMap(ii,jj) + OmegaMap(jj,ii);
            TotalOmega = TotalOmega + OmegaOverlap;
            
            if OmegaOverlap > MaxOmega
                MaxOmega = OmegaOverlap;
                rcMax(1) = ii;
                rcMax(2) = jj;
            end
            if (jj < k) %(k - 1)
                jj = jj + 1;
            else
                ii = ii + 1;
                jj = ii + 1;
            end
            
        end
    end
    
    
    if asympt == 1
        
        while (ii < k)
            
            if fix(ii) == 1
                
                if fix(jj) == 1
                    
                    Di = di(ii,jj,:);
                    Li = li(ii,jj,:);
                    coef = Li - 1.0;
                    ldprod = Li.*Di;
                    const2 = (ldprod.*Di)./coef;
                    ncp = (ldprod./coef).^2;
                    
                    t = sum(const2) + const1(ii,jj);
                    
                    OmegaMap(ii,jj)=ncx2mixtcdf(t, df, coef,ncp,'sigma',sigma,'lim',lim,'tol',tol);
                    
                    % OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
                    
                else
                    
                    OmegaMap(ii,jj) = 0.0;
                    
                end
                
            else
                if fix(jj) == 1
                    
                    OmegaMap(ii,jj) = 0.0;
                    
                else
                    coef = squeeze(li(ii,jj,:)) - 1.0;
                    ncp = zeros(v,1);
                    
                    t = const1(ii,jj);
                    
                    OmegaMap(ii,jj)=ncx2mixtcdf(t, df, coef,ncp,'sigma',sigma,'lim',lim,'tol',tol);
                    
                    % OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
                    
                    
                end
            end
            % rows 508 and 509 of libOverlap.c
            if fix(jj) == 1
                
                if fix(ii) == 1
                    
                    Di = squeeze(di(jj,ii,:));
                    Li = squeeze(li(jj,ii,:));
                    coef = Li - 1.0;
                    ldprod = Li.*Di;
                    const2 = (ldprod.*Di)./coef;
                    ncp = (ldprod./coef).^2;
                    t = sum(const2) + const1(jj,ii);
                    
                    
                    % 	OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
                    
                    OmegaMap(jj,ii)=ncx2mixtcdf(t, df, coef,ncp,'sigma',sigma,'lim',lim,'tol',tol);
                    
                else
                    
                    OmegaMap(jj,ii) = 0.0;
                    
                end
                
            else
                
                if fix(ii) == 1
                    OmegaMap(jj,ii) = 0.0;
                    
                else
                    coef =  squeeze(li(jj,ii,:))- 1.0;
                    ncp = zeros(v,1);
                    
                    t = const1(jj,ii);
                    
                    % OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
                    OmegaMap(jj,ii)=ncx2mixtcdf(t, df, coef,ncp,'sigma',sigma,'lim',lim,'tol',tol);
                    
                end
            end
            
            
            
            
            
            OmegaOverlap = OmegaMap(ii,jj) + OmegaMap(jj,ii);
            TotalOmega = TotalOmega + OmegaOverlap;
            
            if OmegaOverlap > MaxOmega
                MaxOmega = OmegaOverlap;
                rcMax(1) = ii;
                rcMax(2) = jj;
            end
            if (jj < k) %(k - 1)
                jj = jj + 1;
            else
                ii = ii + 1;
                jj = ii + 1;
            end
            
        end
        
        
    end
end
BarOmega=TotalOmega/(0.5*k*(k-1));
end
