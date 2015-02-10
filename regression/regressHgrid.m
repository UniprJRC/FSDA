function [out] = regressHgrid(y,X,sel,varargin)
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
% Last modified 06-Feb-2015


nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

options=struct('intercept',1,'msgiter',0);

if nargin > 3
    
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end
intercept=options.intercept;

%alpha=[1.7217 0.1:0.1:7];
%alpha=3;
%alpha=2;
%alpha=1:0.1:2;
%alpha=1.5;
%alpha=2.9;
%alpha=1;
%gam=[0.000001 0.00001 0.0001 0.001 0.01 0.1:0.1:15]; %  16:1000]; %  6:10000];
%gam=[5.1868e+04 0.001 0.01 0.1  1 10:120 500 1000 5000 10000 50000]; %  6:10000];
alpha=0.1:0.1:4;
gam=[0.001 0.01 0.1 1 10:120 500 1000 5000 10000 50000]; 
%gam=1;
%gam=1;
% const = constant in the scedastc equation
const=0.000001;
const=1;

logL=zeros(length(gam)*length(alpha),3);
%Initialization of gamma
% Z = n-by-r matrix which contains the explanatory variables for
% heteroskedasticity
if size(sel,1)==n
    Z=sel;
else
    % Check if interecept was true
    if intercept==1
        Z=X(:,sel+1);
    else
        Z=X(:,sel);
    end
end

sigma2all=zeros(length(gam)*length(alpha),1);
if intercept==1
betaall=zeros(length(gam)*length(alpha),2);
else
betaall=zeros(length(gam)*length(alpha),1);
end    

ij=1;
for i_alpha=1:length(alpha)
    Zi=real(Z.^alpha(i_alpha));
    for j_gamma=1:length(gam)
        
        sigma2hati=const+Zi*gam(j_gamma);
        
        sqweights = sigma2hati.^(-0.5);
        
        % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
        Xw = bsxfun(@times, X, sqweights);
        yw = y .* sqweights;
        
        % Estimate of beta from (re)weighted regression (RWLS)
        newbeta = Xw\yw;
        
        newres2=(yw-Xw*newbeta).^2;
        % s2 = MLE of sigma2 on transformed data
        sigma2=sum(newres2)/n;
        sigma2all(ij)=sigma2;
        betaall(ij,:)=newbeta';
        % Construct the value of the likelihood
        % loglik=sum(log(sigma2*sigma2hati))+sum(newres2./(sigma2*sigma2hati));
        
        loglik=sum(log(sigma2hati))+n*log(sigma2);
        logL(ij,:)=[alpha(i_alpha) gam(j_gamma) loglik];
        ij=ij+1;
    end
end
[~,indmin]=min(logL(:,3));
% plot(logL(:,3))
out=struct;
out.Beta=betaall(indmin,:);
out.Gamma=logL(indmin,2);
out.alpha=logL(indmin,1);
out.sigma2=sigma2all(indmin);
out.logL=logL;
end
