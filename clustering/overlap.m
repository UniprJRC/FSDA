function [OmegaMap, BarOmega, MaxOmega, rcMax] = overlap(k, p, Pi, Mu, S, tol, lim)
%overlap computes the exact overlap given the parameters of the mixture
%
%  Required input arguments:
%  
%  k  : number of components (not present in R implementation, as it can be
%       derived from the size of Mu)
%  p  : dimensionality (not present in R implementation, as it can be
%       derived from the size of Mu)
%  Pi : mixing proportions
%  Mu : mean vectors (matrix of size k-by-p
%  S  : covariance matrices 3D array of size p-by-p-by-k 
%
%  
%
%  Optional input arguments:
%  
%  tol, lim - parameters for qfc function
%  tol : tolerance (default is 1e-06)
%  lim : scalar = maximum number of integration terms (default is 100000)

% OUTPUT
%
%  - OmegaMap - map of misclassification probabilities
%  - BarOmega - average overlap
%  - MaxOmega - maximum overlap
%  - rcMax - contains the pair of components producing the highest overlap
%
% Copyright 2008-2014. FSDA toolbox
%
%
%<a href="matlab: docsearch('overlap')">Link to the help function</a>
% Last modified 08-Dec-2013
%
% Examples:
%
%
%{
%    Finding exact overlap for the Iris data

    load fisheriris;
    Y=meas;
    Mu=grpstats(Y,species);

    S=zeros(4,4,3);
    S(:,:,1)=cov(Y(1:50,:));
    S(:,:,2)=cov(Y(51:100,:));
    S(:,:,3)=cov(Y(101:150,:));

    pigen=ones(3,1)/3;
    k=3;
    p=4;
    [OmegaMap, BarOmega, MaxOmega, rcMax]=overlap(k,p,pigen,Mu,S)
%}

% Tra paraentisi non capisco come mai l'istruzione che segue non funziona
% grpstats(meas,species,{'cov'})%}

%% Beginning of code

if nargin<5
    error('error: not all input terms have been supplied')
end
if nargin==6
    lim=100000;
end
if nargin ==5
   tol=1e-6;
   lim=100000;
end

[li,di,const1]=ComputePars(p,k,Pi,Mu,S);
fixcl=zeros(k,1);
c=1;
asympt=0;
[OmegaMap, BarOmega, MaxOmega, rcMax]=GetOmegaMap(c, p, k, li, di, const1, fixcl, tol, lim, asympt);

end
