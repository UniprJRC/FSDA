function [out] = regressHart_grid(y,X,sel,varargin)
% grid search to find minimum for ART model
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%
% Last modified 06-Feb-2015

%% Beginning of code

nnargin = nargin;
vvarargin = varargin;
[y,X,n,p] = chkinputR(y,X,nnargin,vvarargin);

alpha=0.1:0.1:4;
gam=[0.001 0.01 0.1 1 10:120 500 1000 5000 10000 50000];

options=struct('intercept',1,'msgiter',0,'const',1,'alpha',alpha,'gam',gam,'interceptSKE',1);

if nargin > 3
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end
% intercept = (boolean scalar) specifying it is necessary to have the
% intercept in the linear regression model
intercept=options.intercept;
% interceptSKE = scalar specifying the value of the intercept in the
% skedastic equation
interceptSKE=options.interceptSKE;

% alpha and gam if supplied by the user provide the grid of values of alpha
% and gamma for which the likelihood must be evaluated
alpha=options.alpha;
gam=options.gam;

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
lalpha=length(alpha);
lgam=length(gam);
lgamlalpha=lgam*lalpha;

sigma2all=zeros(1,lgamlalpha);
if intercept==1
    betaall=zeros(lgamlalpha,2);
else
    betaall=zeros(lgamlalpha,1);
end

% Initialization of Xw;
Xw=X;

ij=1;
sigma2hatall=zeros(n,lgamlalpha);

%logL=zeros(length(gam)*length(alpha),3);

for i_alpha=1:lalpha
    Zialpha=real(Z.^alpha(i_alpha));
    sigma2hat=interceptSKE+Zialpha*gam;
    sigma2hatall(:,(i_alpha-1)*lgam+1:i_alpha*lgam)=sigma2hat;
    % sqweights=sigma2hat.^(-0.5);
    sqweights=exp(-0.5*log(sigma2hat));
    
    for j_gamma=1:lgam
        
        
        % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
        if p==2
            Xw(:,1)=X(:,1).*sqweights(:,j_gamma);
            Xw(:,2)=X(:,2).*sqweights(:,j_gamma);
        elseif p==1
            Xw=X.*sqweights(:,j_gamma);
        else
            Xw = bsxfun(@times, X, sqweights(:,j_gamma));
        end
        yw = y .* sqweights(:,j_gamma);
        
        % Estimate of beta from weighted regression (transformed space)
        beta = Xw\yw;
        
        res2=(yw-Xw*beta).^2;
        % s2 = MLE of sigma2 using transformed data
        sigma2=sum(res2)/n;
        sigma2all(ij)=sigma2;
        betaall(ij,:)=beta';
        
        %         logliktmp=sum(log(sigma2hatall(:,ij)))+n*log(sigma2);
        %         logL(ij,:)=[alpha(i_alpha) gam(j_gamma) logliktmp];
        
        ij=ij+1;
    end
end
loglik=sum(log(sigma2hatall),1)+n*log(sigma2all);
[~,indmin]=min(loglik);

out=struct;
out.Beta=betaall(indmin,:);

% Store values of gamma and alpha which maximized the likelihood
aa=repmat(alpha,lgam,1);
gg=repmat(gam',lalpha,1);
out.Gamma=gg(indmin);
out.alpha=aa(indmin);

% Store estimate of sigma2 (sum of suares of residuals/n in the
% transformed scale)
out.sigma2=sigma2all(indmin);

%Store value of maximized log likelihood
out.logL=loglik(indmin);

% if plots==10
%     aa=reshape(loglik,length(gam),length(alpha));
%     mina=min(min(aa));
%     maxa=max(max(aa));
%     for jj=1:40
%         subplot(5,8,jj)
%         
%         plot(aa(:,jj))
%         xlim([1 lgam])
%         
%         %   plot(gam,aa(:,jj))
%         % xlim([-10000 gam(end)])
%         ylim([mina maxa])
%         title(['\alpha=' num2str(alpha(jj))])
%     end
% end

end