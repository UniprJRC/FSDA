function [OmegaMap, BarOmega, MaxOmega, rcMax]=GetOmegaMap(c, v, k, li, di, const1, fix, tol, lim, asympt, toll)
%GetOmegaMap calculates the map of misclassificaton between groups
%
%<a href="matlab: docsearchFS('GetOmegaMap')">Link to the help function</a>
%
%  Required input arguments:
%
%       c  : inflation parameter for covariance matrices. Scalar. In the
%            case of heterogeneous clusters scalar c is used to correct
%            $li(i,j)$ and $const1(i,j)$. More precisely when hom=0:
%            if $fix(i)=0$, $di(i,j,:)=di(i,j,:)/c^{0.5}$
%                   (because $di(i,j,:)=\Gamma (c \Sigma_i)^{-0.5}(\mu_i-\mu_j)$
%                    where Gamma is the matrix of eigenvectors of
%                    $\Sigma_{j|i}$).
%           if $fix(i)=0$ and $fix(j)=1$,  $li(i,j,:)=c li(i,j,:)$, 
%                    $const1(i,j)=const1(i,j)+v log(c)$
%                   (because the eigenvalues of matrix
%                   $\Sigma_{j|i} = (c^{0.5} \Sigma_i^{0.5}) \Sigma_j^-1
%                   (c^{0.5} \Sigma_i^{0.5})$ are multiplied by $c$.
%                   Similarly $log |c\Sigma_i|= c log(v) + log |\Sigma_i|$.
%           if $fix(i)=0$  and $fix(j)=0$,  $li$ and $const1$ are not changed because
%                   $\Sigma_j|i=(c^0.5 \Sigma_i^{0.5}) (c \Sigma_j) ^-1 (c^0.5
%                   \Sigma_i^{0.5})$, 
%                   $const1(i,j)=log((Pi(j)/Pi(i))^2  c detS(i)/(c  det S(j)))$
%               Data Types - single | double
%       v  : number of variables. Scalar. Dimensionality of the data
%           matrix.
%               Data Types - single | double
%       k  : number of components (groups). Scalar. Scalar associated with
%           the number of groups. 
%               Data Types - single | double
%       li : eigenvalues of matrix $Sji=\Sigma_{j|i}$. 3D array of size
%            k-by-k-by-v. $li(i,j,:)$ is the vector which
%            contains the eigenvalues of matrix $Sji=\Sigma_{j|i}$
%            where $\Sigma_{j|i}  = \Sigma_i^{0.5} \Sigma_j^{-1} \Sigma_i^{0.5}$.
%            Remark: $li(i,j,:)-1$ are the coefficients of the linear
%            combinations of non central $\chi^2$ distribution which is used
%            to compute the probability of overlapping
%            The non centrality parameter  of the non central $\chi^2$ distribution
%            associated with $li(i,j,l)-1$ is $(li(i,j,l)
%            di(i,j,l)/(li(i,j,l)-1))^2$
%               Data Types - single | double
%       di : eigenvector of matrix $\Sigma_{j|i}$. 3D array of size
%            k-by-k-by-v.
%            $di(i,j,l)= \gamma_l' \Sigma_i^{-0.5}(\mu_i-\mu_j)$
%            $l=1, 2, ...,v$, where $\gamma_l$ is the l-th eigenvector coming
%            from the spectral decomposition of matrix $\Sigma_{j|i}$  . Vector
%            di(i,j,:) contains the coefficients of N(0,1) (called W_l in
%            equation 2.2 of Maitra and Melnikov (2010), JCGS)
%               Data Types - single | double
%    const1: Prior and Determinants. Matrix. k-by-k matrix whose element
%            $i,j$ with $i=1,..., k$ and $j=1, ..., k$
%            is equal to $log((Pi(j)/Pi(i))^2 * detS(i)/detS(j))$
%            where $detS(i)$ and $detS(j)$ are, respectively, the
%            determinants of the covariance matrices of groups i and j. See
%            equation (2.2) of Maitra and Melnikov (2010), JCGS.
%            Note that const1(j,i)=-const1(i,j). The elements on the
%            diagonal of matrix const1 are set to 0 (in other words they
%            are not computed).
%           REMARK: li, di and const1 are the parameters needed for computing
%           overlap
%               Data Types - single | double
%      fix : Inflation/Deflation clusters. Vector. Vector of length k containing zeros or ones
%          if fix(j) =1 cluster j does not participate to inflation or
%          deflation. If fix=zeros(k,1) all clusters participate in
%          inflation/deflation.
%          REMARK: this parameter is used just if heterogeneous clusters
%          are used.
%               Data Types - single | double
%      tol : overlap tolerance. Scalar. Error bound for overlap computation default is 1e-06
%               Data Types - single | double
%      lim : maximum number of integration terms. Scalar. default is 1e06.
%               REMARK: Optional parameters tol and lim will be used by
%               function ncx2mixtcdf.m which computes the cdf of a linear
%               combination of non central chi2 r.v.. This is the
%               probability of overlapping
%               Data Types - single | double
%   asympt : flag for regular or asymptotic overlap. Scalar. If asympt ==1 formula
%            of asymptotic overlap is used (see p. 359, paragraph starting
%            with step 3 of Maitra and Melnikov, 2010, JCGS). In this case
%            the misclassification probability $\omega_j|i$ (that is with
%            respect to the centroid and covariance matrix of group i) is
%            given by
%            \[
%            \sum_{l=1}^p (\lambda_l-1) U_l \leq log \pi_j^2 |\Sigma_j| / \pi_i^2 |\Sigma_i|
%            \]
%            where $U_l$ are central chi2 with one degree of freedom.
%            Note that parameter asympt is set to 1 in the preliminary
%            phase just to check whether the average or the maximum overlap
%            can be reached.
%               Data Types - single | double
%
%  Optional input arguments:
%
%     toll : eigenvalues tolerance. Scalar. Tolerance use to declare the eigenvalues of matrix
%           $\Sigma_{j|i}$ equal to 1. The default value is 1e-06
%           Background: the probability of overlapping is the result of two
%           sums. The first is a sum in correspondence of the eigenvalues
%           of matrix $\Sigma_{j|i}$ not equal 1. The second is a sum in
%           correspondence of the eigenvalues of matrix \Sigma_j|i equal 1.
%           Similarly, what is on the right hand side of the probability of
%           overlapping has two sums (for eigenvalues equal or different
%           from 1). Toll specifies when we must consider an eigenvalue
%           equal to 1.
%               Example - 'toll',1e-10
%               Data Types - double
%
%  Output:
%
%    OmegaMap : k-by-k matrix containing map of misclassification
%               probabilities. More precisely, OmegaMap(i,j) is the
%               misclassification probability with respect to cluster i,
%               (that is conditionally on x belonging to cluster i,  which
%               is called $w_{j|i}$) 
%               $(i ~= j)=1, 2, ..., k$.
%    BarOmega : Average overlap. Scalar.
%               BarOmega is computed as (sum(sum(OmegaMap))-k)/(0.5*k(k-1))
%    MaxOmega : Maximum overlap. Scalar. MaxOmega is the
%               maximum of OmegaMap(i,j)+OmegaMap(j,i)
%               $(i ~= j)=1, 2, ..., k$
%       rcMax : Components with highest overlap. Column vector of length
%               equal to 2. It contains the indexes
%               associated with the pair of components producing the
%               highest overlap (largest off diagonal element of matrix
%               OmegaMap).
%
% See also: MixSim, ncx2mixtcdf, restreigen
%
% References:
%
%   Maitra, R. and Melnykov, V. (2010), Simulating data to study performance
%   of finite mixture modeling and clustering algorithms, The Journal of
%   Computational and Graphical Statistics, 2:19, pp. 354-376. (to refer to
%   this publication we will use "MM2010 JCGS")
%
%   Melnykov, V., Chen, W.-C., and Maitra, R. (2012), MixSim: An R Package
%   for Simulating Data to Study Performance of Clustering Algorithms,
%   Journal of Statistical Software, 51:12, pp. 1-25.
%
%   Davies, R. (1980), The distribution of a linear combination of
%   chi-square random variables, Applied Statistics, 29, pp. 323-333.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('GetOmegaMap')">Link to the help function</a>
%
%
%$LastChangedDate::                      $: Date of the last commit
%
%{
    %% Find matrix of misclassification probabilities.
    % Fix number of groups, dimensions and mixing proportions
    k = 4;    % Number of groups
    p = 5;    % Number of dimensions
    Pi= [0.1 0.2 0.4 0.3]; % mixing proportions
    % Mu matrix of means of size k-by-p; (each row is a distinct centroid)
    Mu=randn(k,p);
    % Groups 2 and 3 are far from the other groups
    Mu(2:3,:) = Mu(2:3,:)+10;
    %Mu(3,:) =  Mu(3,:)+2;
    %Mu(4,:) =  Mu(4,:)+3;
    % S= 3D array of dimension p-by-p-by-k containing covariance matrices 
    % of the groups
    S=zeros(p,p,k);
    for j=1:k
        S(:,:,j)=cov(randn(p+10,p));
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

if nargin<11
    toll=1e-6;
end

% tolncx2 = tolerance to use in routine ncx2mixtcdf (which computes cdf of
% linear combinations of non central chi2 distributions)
tolncx2  = tol;

% Omegamap= k-by-k matrix which will contain misclassification probabilities
OmegaMap=eye(k);
rcMax=zeros(2,1);

TotalOmega = 0.0;
MaxOmega = 0.0;

df=ones(v,1);

ii = 1;
jj = 2;

while ii < k
    % check if clusters ii and jj are homogeneous
    hom = 1;
    for kk=1:v
        if abs(li(ii,jj,kk)-li(jj,ii,kk))> eps %1e-14;            
            hom = 0;
            break
        end
    end
    
    if hom == 1   % homogeneous clusters
        if asympt == 0
            Di =  squeeze(di(ii,jj,:)) / sqrt(c);
            Cnst1 = sum(Di.^2);
            
            % t is the value in which the cdf of the mixture of
            % non central chi2 must be evaluated.
            t = const1(ii,jj) - Cnst1;
            sigma = 2 * sqrt(Cnst1);
            
            % This is nothing but the cdf of the standard normal of
            % \Phi( \frac{log \pi_j^2/\pi_i^2 - \sum_{l=1}^v \delta_l^2}
            %            { 2 \sqrt{ \sum_{l=1}^v \delta_l^2 }=
            % \Phi( \frac{log \pi_j/\pi_i}{\sqrt{ \sum_{l=1}^v \delta_l^2}
            %       -0.5 sqrt{\sum_{l=1}^v \delta_l^2} =
            % Using the symbols in this code this means we simply have to
            % compute \Phi( t / sigma)
            % See also Melnikov, Chen and Maitra, JSS 2012, p. 4
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
        else % asympt==1
            % The value of asymptotic overlap for homogeneous clusters
            % omega_{ij}^\infty = omega_{j|i}^\infty+omega_{i|j}^\infty =1 for
            % any mixing proportions \pi_i and \pi_j
            % See for example Melnikov, Chen and Maitra (2012) , JSS p. 4
            
            
            % if \pi_j > \pi_i (const1(ii,jj) > 0)
            % \omega_{j|i}= OmegaMap(ii,jj) =1
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
            
        end
    else % non homogeneous clusters
        sigma = 0.0;
        if asympt == 0
            
            lij=squeeze(li(ii,jj,:));
            dij=squeeze(di(ii,jj,:));
            
            lijne1=abs(lij-1)>toll;
            
            
            % if sum(lijne1)<v there are eigenvalues which are  = 1
            if sum(lijne1)<v
                eigeq1=1;
                lijeq1= (~lijne1);
            else
                eigeq1=0;
            end
            
            if fix(ii) == 1
                % Di = squeeze(di(ii,jj,lijne1));
                Di = dij(lijne1);
                if eigeq1==1
                    Dieq1 = dij(lijeq1);
                end
                
                if fix(jj) == 1
                    % Li = squeeze(li(ii,jj,:));
                    Li = lij(lijne1);
                    Cnst1 = const1(ii,jj);
                else
                    % Li = squeeze(li(ii,jj,:))/ c;
                    Li = lij(lijne1) / c;
                    Cnst1 = const1(ii,jj) - v * log(c);
                end
                
            else
                
                
                % Di = squeeze(di(ii,jj,:)) / sqrt(c);
                Di = dij(lijne1) / sqrt(c);
                if eigeq1==1
                    Dieq1 = dij(lijeq1)/sqrt(c);
                end
                
                if fix(jj) == 1
                    
                    % Li = c * squeeze(li(ii,jj,:));
                    Li = c * lij(lijne1);
                    Cnst1 = const1(ii,jj) + v * log(c);
                else
                    % Li = squeeze(li(ii,jj,:));
                    Li = lij(lijne1);
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
            
            if eigeq1==1
                sDieq1sq=sum(Dieq1.^2);
                t= t - sDieq1sq;
                % sigmaall= sigma+ 2*sum(Dieq1);
                sigmaall= sigma+ 2*sqrt(sDieq1sq);
            else
                sigmaall=sigma;
            end
            
            OmegaMap(ii,jj)=ncx2mixtcdf(t, df(1:sum(lijne1)), coef, ncp,'sigma',sigmaall,'lim',lim,'tol',tolncx2);
            
            lji=squeeze(li(jj,ii,:));
            dji=squeeze(di(jj,ii,:));
            ljine1=abs(lji-1)>toll;
            
            if sum(ljine1)<v
                eigeq1=1;
                ljieq1= (~ljine1);
            else
                eigeq1=0;
            end
            
            if fix(jj) == 1
                % Di = squeeze(di(jj,ii,:));
                Di = dji(ljine1);
                if eigeq1==1
                    Dieq1 = dji(ljieq1);
                end
                
                if fix(ii) == 1
                    % Li = squeeze(li(jj,ii,:));
                    Li =  lji(ljine1);
                    Cnst1 = const1(jj,ii);
                else
                    % Li = squeeze(li(jj,ii,:)) / c;
                    Li =  lji(ljine1)/c;
                    Cnst1 = const1(jj,ii) - v * log(c);
                end
                
            else
                
                % Di = squeeze(di(jj,ii,:)) / sqrt(c);
                Di = dji(ljine1) / sqrt(c);
                if eigeq1==1
                    Dieq1 = dji(ljieq1) / sqrt(c);
                end
                
                if fix(ii) == 1
                    % Li = c * squeeze(li(jj,ii,:));
                    Li = c * lji(ljine1);
                    
                    Cnst1 = const1(jj,ii) + v * log(c);
                else
                    % Li = squeeze(li(jj,ii,:));
                    Li = lji(ljine1);
                    
                    Cnst1 = const1(jj,ii);
                end
                
            end
            
            coef = Li - 1.0;
            ldprod = Li.*Di;
            const2 = (ldprod.*Di)./ coef;
            ncp = (ldprod./coef).^2;
            t = sum(const2)+ Cnst1;
            
            if eigeq1==1
                sDieq1sq=sum(Dieq1.^2);
                t = t - sDieq1sq;
                % sigmaall= sigma+ 2*sum(Dieq1);
                sigmaall= sigma+ 2*sqrt(sDieq1sq);
            else
                sigmaall=sigma;
            end
            OmegaMap(jj,ii)=ncx2mixtcdf(t, df(1:sum(ljine1)), coef, ncp,'sigma',sigmaall,'lim',lim,'tol',tolncx2);
            
            OmegaOverlap = OmegaMap(ii,jj) + OmegaMap(jj,ii);
            TotalOmega = TotalOmega + OmegaOverlap;
            
            if OmegaOverlap > MaxOmega
                MaxOmega = OmegaOverlap;
                rcMax(1) = ii;
                rcMax(2) = jj;
            end
        else % asympt =1
            lij=squeeze(li(ii,jj,:));
            dij=squeeze(di(ii,jj,:));
            
            
            lijne1=abs(lij-1)>toll;
            
            % if sum(lijne1)<v there are eigenvalues which are  = 1
            if sum(lijne1)<v
                eigeq1=1;
                lijeq1= (~lijne1);
            else
                eigeq1=0;
            end
            
            
            if fix(ii) == 1
                
                if fix(jj) == 1
                    
                    % Di = di(ii,jj,:);
                    Di = dij(lijne1);
                    
                    if eigeq1==1
                        Dieq1 = dij(lijeq1);
                    end
                    %Li = li(ii,jj,:);
                    Li = lij(lijne1);
                    
                    
                    coef = Li - 1.0;
                    ldprod = Li.*Di;
                    const2 = (ldprod.*Di)./coef;
                    ncp = (ldprod./coef).^2;
                    
                    t = sum(const2) + const1(ii,jj);
                    
                    if eigeq1==1
                        sDieq1sq=sum(Dieq1.^2);
                        t = t - sDieq1sq;
                        % sigmaall= sigma+ 2*sum(Dieq1);
                        sigmaall= sigma+ 2*sqrt(sDieq1sq);
                    else
                        sigmaall=sigma;
                    end
                    OmegaMap(ii,jj)=ncx2mixtcdf(t, df(1:sum(lijne1)), coef, ncp,'sigma',sigmaall,'lim',lim,'tol',tolncx2);
                    
                    % OmegaMap[i][j] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
                    
                else
                    % Because in this case Const1 is - \infty
                    OmegaMap(ii,jj) = 0.0;
                    
                end
                
            else
                if fix(jj) == 1
                    % Because in this case
                    % Pr(\infty * chi^2_1 (central) < (Const1 + log(\infty))=
                    % Pr(chi^2_1 (central) < 0)=0
                    OmegaMap(ii,jj) = 0.0;
                    
                else
                    % If c \rightarrow \infty (asympt =1) and cluster are
                    % heterogeneous (hom=0) the probability of overlapping
                    % is \sum_{l=1}^k li(ii,jj,l) U_l \leq const1(ii,jj) = t
                    % where U_l are independent (central) \chi^2_1 (with
                    % one degree of freedom)
                    Li = lij(lijne1);
                    
                    % coef = squeeze(li(ii,jj,:)) - 1.0;
                    coef = Li - 1.0;
                    
                    ncp = zeros(sum(lijne1),1);
                    
                    t = const1(ii,jj);
                    
                        OmegaMap(ii,jj)=ncx2mixtcdf(t, df(1:sum(lijne1)), coef, ncp,'sigma',sigma,'lim',lim,'tol',tolncx2);
                                
                end
            end
            
            lji=squeeze(li(jj,ii,:));
            dji=squeeze(di(jj,ii,:));
            ljine1=abs(lji-1)>toll;
            
            if sum(ljine1)<v
                eigeq1=1;
                ljieq1= (~ljine1);
            else
                eigeq1=0;
            end
            
            
            % rows 508 and 509 of libOverlap.c
            if fix(jj) == 1
                
                if fix(ii) == 1
                    
                    %  Di = squeeze(di(jj,ii,:));
                    Di = dji(ljine1);
                    
                    if eigeq1==1
                        Dieq1 = dji(ljieq1);
                    end
                    %Li = squeeze(li(jj,ii,:));
                    Li = lji(ljine1);
                    
                    coef = Li - 1.0;
                    ldprod = Li.*Di;
                    const2 = (ldprod.*Di)./coef;
                    ncp = (ldprod./coef).^2;
                    t = sum(const2) + const1(jj,ii);
                    
                    if eigeq1==1
                        sDieq1sq=sum(Dieq1.^2);
                        t = t - sDieq1sq;
                        % sigmaall= sigma+ 2*sum(Dieq1);
                        sigmaall= sigma+ 2*sqrt(sDieq1sq);
                        
                    else
                        sigmaall=sigma;
                    end
                    % 	OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
                    OmegaMap(jj,ii)=ncx2mixtcdf(t, df(1:sum(ljine1)), coef,ncp,'sigma',sigmaall,'lim',lim,'tol',tolncx2);
                    
                else
                    
                    OmegaMap(jj,ii) = 0.0;
                    
                end
                
            else
                
                if fix(ii) == 1
                    OmegaMap(jj,ii) = 0.0;
                    
                else
                    % coef =  squeeze(li(jj,ii,:))- 1.0;
                    Li = lji(ljine1);
                    coef = Li - 1.0;
                    ncp = zeros(sum(ljine1),1);
                    t = const1(jj,ii);
                    % OmegaMap[j][i] = qfc(coef, ncp, df, &p, &sigma, &t, &lim, &acc, trace, &ifault);
                    OmegaMap(jj,ii)=ncx2mixtcdf(t, df(1:sum(ljine1)), coef, ncp,'sigma',sigma,'lim',lim,'tol',tolncx2);
                    
                end
            end
            
            
            OmegaOverlap = OmegaMap(ii,jj) + OmegaMap(jj,ii);
            TotalOmega = TotalOmega + OmegaOverlap;
            
            if OmegaOverlap > MaxOmega
                MaxOmega = OmegaOverlap;
                rcMax(1) = ii;
                rcMax(2) = jj;
            end
        end
    end
    if jj < k
        jj = jj + 1;
    else
        ii = ii + 1;
        jj = ii + 1;
    end
    
end

BarOmega=TotalOmega/(0.5*k*(k-1));
end