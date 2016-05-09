function [OmegaMap, BarOmega, MaxOmega, rcMax]=GetOmegaMapv1(c, Pi, Mu, S, fix, asympt)
%GetOmegaMapv1 calculates the map of misclassificaton between groups when v=1
%
%  Required input arguments:
%
%       c  : scalar, inflation parameter for covariance matrices. In the
%            case of heterogeneous clusters scalar c is used to correct
%            the elements of S for which the corresponding element of
%            vector fix is 0
%      Pi  : vector of length k containing the mixing proportions.
%            Clearly, sum(Pi)=1.
%      Mu  : vector of length k containing the k centers of the groups
%       S  : vector of length k containing the variances of the k groups
%      fix : vector of length k containing zeros or ones
%            if fix(j) =1 cluster j does not participate to inflation or
%            deflation. If fix=zeros(k,1) all clusters participate in
%            inflation/deflation
%            This parameter is used just if heterogeneous clusters are used
%            (that is when max(S)>min(S))
%   asympt : flag for regular or asymptotic overlap. If asypmt ==1 formula
%            of asymptotic overlap is used. In this case
%            the misclassification probability \omega_j|i (that is with
%            respect to the mean and variance of group i) is
%            computed assuming that the variance of group i goes to infinity
%            Note that parameter asympt is set to 1 in the preliminary
%            phase just to check whether the average or the maximum overlap
%            can be reached.
%
%
%  Output:
%
%    OmegaMap : k-by-k matrix containing map of misclassification
%               probabilities. More precisely, OmegaMap(i,j) is the
%               misclassification probability with respect to cluster i,
%               (that is conditionally on x belonging to cluster i,  which
%               is called w_{j|i})
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
% Copyright 2008-2015.
% Written by FSDA team
%
% Last modified 06-Feb-2015
%{
    % COMPARE OLD AND NEW FUNCTION
    rng(65123,'twister');
    k=4; % Number of groups
    p=1; % Number of dimensions
    Pi=[0.1 0.2 0.4 0.3]; % mixing proportions
    % Mu = vector of length k contains the means of the k groups
    Mu=randn(k,1);
    % Groups 2 and 3  is far from the other groups
    Mu(2:3)=Mu(2:3)+2.2;
    % S= 3D array of dimension 1-by-1-by-k containing the variances of
    % the k groups
    S=zeros(p,p,k);
    % Generate groups with unequal variance
    hom=0;
    if hom ==1
    S(1,1,:)=abs(randn(1,1));
    else
    S(1,1,:)=abs(randn(k,1));
    end
    asympt = 0;
    % Inflate variances using parameter c
    c = 2.3;
    % Cluster 2 and 3 do not participate to the inflation/deflation process with
    % parameter c
    fixcl=zeros(k,1);
    fixcl(2)=1;
    fixcl(3)=1;

    % Use traditional functions
    [li, di, const1]=ComputePars(p, k, Pi, Mu, S);
    tol=1e-8;
    lim=1e+07;
    [OmegaMap1,  BarOmega, MaxOmega, rcMax] = ...
        GetOmegaMap(c, p, k, li, di, const1, fixcl, tol, lim, asympt);

    % New function
    S=squeeze(S);
    [OmegaMapv1, BarOmegav1, MaxOmegav1, rcMaxv1]=GetOmegaMapv1(c, Pi, Mu, S, fixcl, asympt);
    disp('OmegaMap using functions ComputePars and GetOmegaMap setting v=1')
    disp(OmegaMap1)
    disp('OmegaMap using new ad hoc function GetOmegaMapv1')
    disp(OmegaMapv1)
%}

%% Beginning of code

verMatlab=verLessThan('matlab','8.2.0');
% newncx2cdf is  a boolean. If it is true new version ncx2cdf is
% used to compute the probability of the upper tail, else if it is false
% old function ncx2cdf is used with argument upper.
if verMatlab ==1
    newncx2cdf=false;
else
    newncx2cdf=true;
end

k=length(Mu);
% Omegamap= k-by-k matrix which will contain misclassification probabilities
OmegaMap=eye(k);
% rcMax vector of length 2 which will contains the indexes of the two
% groups which show the highest overlap
rcMax=zeros(2,1);

TotalOmega = 0.0;
MaxOmega = 0.0;


% Initialize ii and jj
ii = 1;
jj = 2;

S1=S;

if asympt ==0
    % Parameter fix is just used in presence of heterogeneous clusters
    % THIS MUST BE DISCUSSED
    % THIS MUST BE DISCUSSED
    % THIS MUST BE DISCUSSED
    % THIS MUST BE DISCUSSED
    % THIS MUST BE DISCUSSED
    
    if max(S)>min(S)
        S1(fix==0)=c*S1(fix==0);
    else
        S1=c*S1;
    end
end


while ii < k
    s2j=S1(jj);
    s2i=S1(ii);
    pii=Pi(ii);
    pij=Pi(jj);
    
    mui=Mu(ii);
    muj=Mu(jj);
    xij=s2j*(mui-muj)^2/(s2i-s2j)^2 + log((pij/pii)^2*(s2i/s2j))*s2j/(s2i-s2j);
    deltaij=s2i*(((mui-muj)/(s2i-s2j))^2);
    
    xji=s2i*(mui-muj)^2/(s2i-s2j)^2 + log((pii/pij)^2*(s2j/s2i))*s2i/(s2j-s2i);
    deltaji=s2j*(((mui-muj)/(s2i-s2j))^2);
    
    if asympt == 0
        
        if s2i>s2j
            OmegaMap(ii,jj) = ncx2cdf(xij,1,deltaij);
            if newncx2cdf
                OmegaMap(jj,ii) = ncx2cdf(xji,1,deltaji,'upper');
            else
                OmegaMap(jj,ii) = 1- ncx2cdf(xji,1,deltaji);
            end
        elseif s2i == s2j
            
            OmegaMap(ii,jj)= normcdf(-0.5*abs(mui-muj)/sqrt(s2i)+sqrt(s2i)*log(pij/pii)/abs(mui-muj));
            OmegaMap(jj,ii)= normcdf(-0.5*abs(mui-muj)/sqrt(s2i)+sqrt(s2i)*log(pii/pij)/abs(mui-muj));
        else
            if newncx2cdf
                OmegaMap(ii,jj) = ncx2cdf(xij,1,deltaij,'upper');
            else
                OmegaMap(ii,jj) = 1-ncx2cdf(xij,1,deltaij);
            end
            OmegaMap(jj,ii) = ncx2cdf(xji,1,deltaji);
        end
    else % use formulae for asymptotic overlap
        
        if  s2i == s2j % Asymptotic overlap with equal variances
            const1=s2j/s2i;
            if const1 > 0
                OmegaMap(ii,jj) = 1;
                OmegaMap(jj,ii) = 0;
            end
            
            if const1 < 0
                OmegaMap(ii,jj) = 0;
                OmegaMap(jj,ii) = 1;
            end
            
            if const1 == 0
                OmegaMap(ii,jj) = 0.5;
                OmegaMap(jj,ii) = 0.5;
            end
        else % Asymptotic overlap with unequal variances
            if fix(ii)==1 && fix(jj) ==1
                % If both groups do not participate to the inflation
                % deflation process than use the usual formulae
                if s2i>s2j
                    OmegaMap(ii,jj) = ncx2cdf(xij,1,deltaij);
                    if newncx2cdf
                        OmegaMap(jj,ii)=ncx2cdf(xji,1,deltaji,'upper');
                    else
                        OmegaMap(jj,ii)=1-ncx2cdf(xji,1,deltaji);
                    end
                elseif s2i == s2j
                    
                    OmegaMap(ii,jj)=normcdf(-0.5*abs(mui-muj)/sqrt(s2i)+sqrt(s2i)*log(pij/pii)/abs(mui-muj));
                    OmegaMap(jj,ii)=normcdf(-0.5*abs(mui-muj)/sqrt(s2i)+sqrt(s2i)*log(pii/pij)/abs(mui-muj));
                else
                    if newncx2cdf
                        OmegaMap(ii,jj) = ncx2cdf(xij,1,deltaij,'upper');
                    else
                        OmegaMap(ii,jj) = 1 - ncx2cdf(xij,1,deltaij);
                    end
                    OmegaMap(jj,ii)=ncx2cdf(xji,1,deltaji);
                end
                
            elseif fix(ii)+fix(jj) ==1 % just the variance of one group goes to infinity
                OmegaMap(ii,jj) = 0;
                OmegaMap(jj,ii) = 0;
            else  % fixcl(ii)==0 && fixcl(jj) ==0
                % central chi2 is used
                if s2i>s2j
                    OmegaMap(ii,jj)= chi2cdf((1/(s2i/s2j-1))*log((pij/pii)^2*s2i/s2j),1);
                    OmegaMap(jj,ii)= chi2cdf((1/(s2j/s2i-1))*log((pii/pij)^2*s2j/s2i),1,'upper');
                else
                    OmegaMap(ii,jj)= chi2cdf((1/(s2i/s2j-1))*log((pij/pii)^2*s2i/s2j),1,'upper');
                    OmegaMap(jj,ii)= chi2cdf((1/(s2j/s2i-1))*log((pii/pij)^2*s2j/s2i),1);
                end
            end
            
        end
        
    end
    
    
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

BarOmega=TotalOmega/(0.5*k*(k-1));
end
%FScategory:CLUS-MixSim