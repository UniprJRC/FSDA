function S=genSigmaGPCM(v, k, modeltype)
%genSigmaGPCM generates covariance matrix for the 14 Gaussian Parsimonious Clustering Models
%
%<a href="matlab: docsearchFS('genSigmaGPCM')">Link to the help function</a>
%
%
%
%  Required input arguments:
%
%            v: number of dimensions (variables). Scalar.
%               Desired number of variables.
%               Data Types - int16|int32|int64|single|double
%
%            k: number of groups (components). Scalar.
%               Desired number of groups.
%               Data Types - int16|int32|int64|single|double
%
%    modeltype: type of Gaussian Parsimonious Clustering Model. Character.
%               A 3 letter word in the set:
%               'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',
%               'EVI','VEI','EEI','VII','EII'
%               Data Types - Character
%
%
%  Optional input arguments:
%
%
%  Output:
%
%
%             S  : v-by-v-by-k array containing covariances for the k
%                  groups.
%
%
%
%  More About:
%
%  Generate covariance matrices from the 14 parsimonious Gaussian
%  clustering models (GPCM).
%
%
% See also: MixSim, restreigen, restrdeter
%
% References:
%
% Celeux, G., Govaert, G. (1995), Gaussian parsimonious clustering models,
% "Pattern Recognition", 28, pp. 781-793.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('genSigmaGPCM')">Link to the help function</a>
%
%$LastChangedDate:: 2018-09-15 00:27:12 #$: Date of the last commit
%

% Examples:
%
%{
    %%  Covariance matrices contours for the 14 models.
    % Two dimensions
    v=2;
    % 3 groups
    k=3;

    models={'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',...
        'EVI','VEI','EEI','VII','EII'};
    % Specify the colors for the ellipses
    col='rbk';
    % if withseed is true the same plot is always obtained otherwise every time
    % a different plot is obtained
    withseed=true;
    close all
    % These numbers are those which better exemplify the caractheristics of the
    % 14 specifications.
    seeds=[100 20 12 209 51 6 76 8 9 29 111 12 130 14];

    for j=1:length(models)
        if withseed==true
            rng(seeds(j))
        end
        modeltype=models{j};
        S=genSigmaGPCM(v, k, modeltype);
        subplot(4,4,j)
        hold('on')
        for i=1:k
            ellipse(zeros(v,1), S(:,:,i),0.95,col(i));
        end
        axis equal
        legend('off')
        title(modeltype)
    end
%}

%{
    %%  Covariance matrices contours for the 14 models (with different locations).
    % Two dimensions
    v=2;
    % 3 groups
    k=3;

    models={'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',...
        'EVI','VEI','EEI','VII','EII'};
    % Specify the colors for the ellipses
    col='rbk';
    % if withseed is true the same plot is always obtained otherwise every time
    % a different plot is obtained
    withseed=true;
    close all
    % These numbers are those which better exemplify the caractheristics of the
    % 14 specifications.
    seeds=[100 20 12 209 51 6 76 8 9 29 111 12 130 14];

    for j=1:length(models)
        if withseed==true
            rng(seeds(j))
        end
        modeltype=models{j};
        S=genSigmaGPCM(v, k, modeltype);
        subplot(4,4,j)
        hold('on')
        for i=1:k
            cen=zeros(v,1)+i*2;
            ellipse(cen, S(:,:,i),0.95,col(i));
        end
        axis equal
        legend('off')
        title(modeltype)
    end
%}

%% Beginning of code
% S = 3d array which contains the covariance matriced of the groups
S=zeros(v,v,k);
coef=1;
% Determinant
lmd=coef*rand(k,1);
% 3D array which will contain shape matrices
GAM3D=zeros(v,v,k);
% 3D array which will contain rotation matrices
OMG3D=GAM3D;

n=10;
rho=2*rand(1)-1;
chrho=chol(gallery('kms',v,rho));

for j=1:k
    gam=coef*rand(v,1);
    gam=gam./(prod(gam)^(1/v));
    GAM3D(:,:,j)=diag(gam);
    
    u=coef*randn(n,v)*chrho;
    covu=cov(u);
    
    [OMG,~] = eig(covu);
    OMG3D(:,:,j)=OMG;
end


if strcmp(modeltype,'VVE')
    
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Equal orientation
        OMG3Dj=OMG3D(:,:,1);
        S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'EVE')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Equal orientation
        OMG3Dj=OMG3D(:,:,1);
        S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'VVV')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Variable orientation
        OMG3Dj=OMG3D(:,:,j);
        S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'EVV')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Variable orientation
        OMG3Dj=OMG3D(:,:,j);
        S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'VEE')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Equal orientation
        OMG3Dj=OMG3D(:,:,1);
        S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'EEE')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Equal orientation
        OMG3Dj=OMG3D(:,:,1);
        S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'VEV')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Variable orientation
        OMG3Dj=OMG3D(:,:,j);
        S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'EEV')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Variable orientation
        OMG3Dj=OMG3D(:,:,j);
        S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'VVI')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Axis-Aligned
        S(:,:,j)=(lmdj^(1/v)) *   GAM3Dj ;
    end
    
elseif    strcmp(modeltype,'EVI')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Axis-Aligned
        S(:,:,j)=(lmdj^(1/v)) *   GAM3Dj ;
    end
    
elseif    strcmp(modeltype,'VEI')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Axis-Aligned
        S(:,:,j)=(lmdj^(1/v)) *   GAM3Dj ;
    end
    
elseif    strcmp(modeltype,'EEI')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Axis-Aligned
        S(:,:,j)=(lmdj^(1/v)) *   GAM3Dj ;
    end
    
elseif    strcmp(modeltype,'VII')
    eyev=eye(v);
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Spherical shape
        S(:,:,j)=(lmdj^(1/v)) *   eyev ;
    end
    
elseif    strcmp(modeltype,'EII')
    eyev=eye(v);
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Spherical shape
        S(:,:,j)=(lmdj^(1/v)) *   eyev ;
    end
else
    error('FSDA:genSigmaGPCM:WrongModelSupplied','Model supplied id not recognized');
end

end
%FScategory:CLUS-RobClaMULT