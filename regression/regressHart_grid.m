function [out] = regressHart_grid(y,X,Z,varargin)
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
theta=[0.001 0.01 0.1 1 1.71 10:120 500 1000 5000 10000 50000];

options=struct('intercept',1,'msgiter',0,'const',1,...
    'alpha',alpha,'theta',theta,'plots',0);

if nargin > 3
    % Write in structure 'options' the options chosen by the user
    for i=1:2:length(varargin);
        options.(varargin{i})=varargin{i+1};
    end
end
% intercept = (boolean scalar) specifying it is necessary to have the
% intercept in the linear regression model
intercept=options.intercept;

% alpha and theta if supplied by the user provide the grid of values of alpha
% and gamma for which the likelihood must be evaluated
alpha=options.alpha;
theta=options.theta;

% Z = n-by-r matrix which contains the explanatory variables for
% heteroskedasticity
if size(Z,1)~=n
    % Check if interecept was true
    if intercept==1
        Z=(X(:,Z+1));
    else
        Z=(X(:,Z));
    end
end
lalpha=length(alpha);
ltheta=length(theta);
lthetalalpha=ltheta*lalpha;

sigma2all=zeros(1,lthetalalpha);
betaall=zeros(p,lthetalalpha);
    
% Initialization of Xw;
Xw=X;

ij=1;
omegahatall=zeros(n,lthetalalpha);

% logLtmp=zeros(ltheta*lalpha,3);

for i_alpha=1:lalpha
    Zialpha=real(Z.^alpha(i_alpha));
    omegahat=1+Zialpha*theta;
    omegahatall(:,(i_alpha-1)*ltheta+1:i_alpha*ltheta)=omegahat;
    % sqweights=omegahat.^(-0.5);
    sqweights=exp(-0.5*log(omegahat));
    
    for j_theta=1:ltheta
        
        
        % Xw = [X(:,1) .* sqweights X(:,2) .* sqweights ... X(:,end) .* sqweights]
        if p==2
            Xw(:,1)=X(:,1).*sqweights(:,j_theta);
            Xw(:,2)=X(:,2).*sqweights(:,j_theta);
        elseif p==1
            Xw=X.*sqweights(:,j_theta);
        else
            Xw = bsxfun(@times, X, sqweights(:,j_theta));
        end
        yw = y .* sqweights(:,j_theta);
        
        % Estimate of beta from weighted regression (transformed space)
        beta = Xw\yw;
        
        res2=(yw-Xw*beta).^2;
        % sigma2 = MLE of sigma2 using transformed data
        sigma2=sum(res2)/n;
        sigma2all(ij)=sigma2;
        betaall(:,ij)=beta;
                
               %  logliktmp=sum(log(omegahat(:,j_theta)))+n*log(sigma2);
               %  logLtmp(ij,:)=[alpha(i_alpha) theta(j_theta) logliktmp];
        
        ij=ij+1;
    end
end
loglik=sum(log(omegahatall),1)+n*log(sigma2all);
[~,indmin]=min(loglik);

% [~,indmintmp]=min(logLtmp(:,3));


out=struct;
% Store in a column vector of estimate regression coefficients 
out.Beta=betaall(:,indmin);

% Store values of gamma and alpha which maximized the likelihood
aa=repmat(alpha,ltheta,1);
gg=repmat(theta',lalpha,1);
out.Gamma= [log(gg(indmin)); aa(indmin)];

% Store estimate of sigma2 (sum of squares of residuals/n in the
% transformed scale)
out.sigma2=sigma2all(indmin);

out.GammaOLD=gg(indmin);
out.alphaOLD=aa(indmin);

%Store value of maximized log likelihood
out.logL=loglik(indmin);

plots=options.plots;

if plots==1
    aa=reshape(loglik,length(theta),length(alpha));
    mina=min(min(aa));
    maxa=max(max(aa));
    for jj=1:40
        subplot(5,8,jj)
        
        plot(aa(:,jj))
        xlim([1 ltheta])
        
        %   plot(gam,aa(:,jj))
        % xlim([-10000 gam(end)])
        ylim([mina maxa])
        title(['\alpha=' num2str(alpha(jj))])
    end
end

end