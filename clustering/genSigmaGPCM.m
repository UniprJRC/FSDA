function S=genSigmaGPCM(v, k, pa)
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
%      pa  : Constraints to apply and model specification. Structure.
%            Structure containing the following fields:
%             pa.pars= type of Gaussian Parsimonious Clustering Model. Character.
%               A 3 letter word in the set:
%               'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',
%               'EVI','VEI','EEI','VII','EII'.
%               The field pa.pars is compulsory. All the other fields are
%               non necessary. If they are not present they are set to
%               their default values.
%             pa.exactrestriction = boolean. If pa.exactrestriction is true
%               the covariance matrices have to be generated with the exact
%               values of the restrictions specified in pa.cdet, pa.shw and
%               pa.swb. In order to reach this purpose, this procedure
%               ensures to generate shape matrices in which at least one
%               ratio of the elements inside each component is greater or
%               equal than that specified is pa.shw and at least one ratio
%               among the ordered elements of each shape matrix is greater or
%               equal tahn pa.shb. The successive application of routine
%               restrSgimaGPCM guarrantes that the inequalites become
%               equalities. The default value of pa.exaxtrestriction is
%               false therefore covariance matrices are generated without
%               implying any constraint.  If pa.exactrestriction is true
%               covariance matrices are generated with:
%               1) max ratio between determinants equal to pa.cdet (if
%               pa.cdet is specififed);
%               2) at least a  ratio between the elements of each shape matrix greater or equal than
%               pa.shw (if pa.shw is specififed);
%               3) max ratio between the ordered elements of each shape
%               matrices greater or equal than pa.shb).
%             pa.cdet = scalar in the interval [1 Inf) which specifies the
%               the restriction which has to be applied to the determinants.
%               This field is used just if pa.exactrestriction is true.
%             pa.shw = scalar in the interval [1 Inf) which specifies the
%               the restriction which has to be applied to the elements of
%               the shape matrices inside each group.  This field is used
%               just if pa.exactrestriction is true.
%             pa.shb = scalar in the interval [1 Inf) which specifies the
%               the restriction which has to be applied to the ordered elements of
%               the shape matrices across groups. This field is used
%               just if pa.exactrestriction is true.
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
% Copyright 2008-2021.
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
    pa=struct;
    models={'VVE','EVE','VVV','EVV','VEE','EEE','VEV','EEV','VVI',...
        'EVI','VEI','EEI','VII','EII'};
    % Specify the colors for the ellipses
    col='rbk';
    % if withseed is true the same plot is always obtained otherwise every time
    % a different plot is obtained
    withseed=true;
    close all
    % These numbers are those which better exemplify the characteristics of the
    % 14 specifications.
    seeds=[100 20 12 209 51 6 76 8 9 29 111 12 130 14];
    pa=struct;
    for j=1:length(models)
        if withseed==true
            rng(seeds(j))
        end
        pa.pars=models{j};
        S=genSigmaGPCM(v, k, pa);
        subplot(4,4,j)
        hold('on')
        for i=1:k
            ellipse(zeros(v,1), S(:,:,i),0.95,col(i));
        end
        axis equal
        legend('off')
        title(pa.pars)
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
    pa=struct;
    for j=1:length(models)
        if withseed==true
            rng(seeds(j))
        end
        modeltype=models{j};
        pa.pars=modeltype;
        S=genSigmaGPCM(v, k, pa);
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

exactrestr=isfield(pa,'exactrestriction')  && pa.exactrestriction==true;

% S = 3d array which contains the covariance matriced of the groups
S=zeros(v,v,k);
coef=1;
% Determinant
lmd=coef*rand(k,1);

if exactrestr==true
    onedivv=1/v;
    if isfield(pa,'cdet')
        cdet=pa.cdet;
    else
        error('FSDA:genSigmaGPCM:WrongModelSupplied',['Exact restriction ' ...
            'for cdet has been chosen but cdet has '...
            'not been supplied inside input structure pa']);
    end
    
    
    condcdet=false;
    itercondcdet=0;
    maxitercondcdet=100000;
    while condcdet==false && itercondcdet<maxitercondcdet
        itercondcdet=itercondcdet+1;
        if max(lmd)/min(lmd)>=cdet^(1/v)
            condcdet=true;
        else
            
            % Generate beta random numbers with parameters 0.5 and
            % 0.5 in order to increase the probability of large
            % values for the ratio max/min.
            % Generate gamma random values and take ratio of the
            % first to the sum.
            g1 = randg(0.5,k,1);
            g2 = randg(0.5,k,1);
            lmd = g1 ./ (g1 + g2);
            % The 3 above instructions are much faster than
            % lmd=betarnd(0.5,0.5,v,1);
            
            % gam=sort(rand(v,1));
        end
    end
    if itercondcdet==maxitercondcdet
        error('FSDA:restrshapeExact:Wwrongconstraints','cdet is too large please decrease cdet');
    else
        % lmd=(lmd);
    end
    niini=ones(k,1);
    lmd = restreigen(lmd',niini, cdet^(1/v),1e-12,2);    %TODOTODO
end


% 3D array which will contain shape matrices
GAM3D=zeros(v,v,k);
% 3D array which will contain rotation matrices
OMG3D=GAM3D;
GAM2D=zeros(v,k);

n=10;
rho=2*rand(1)-1;
chrho=chol(gallery('kms',v,rho));


if exactrestr==true
    if isfield(pa,'shb') && isfield(pa,'shw')
        shw=pa.shw;
        shb=pa.shb;
    else
        error('FSDA:genSigmaGPCM:WrongModelSupplied',['Exact restriction ' ...
            'for shw and shb have been chosen but shb and shw have '...
            'not been supplied inside input structure pa']);
    end
    
    Sh=restrshapeExact(k, v, shw, shb);
    for j=1:k
        GAM3D(:,:,j)=diag(Sh(:,j));
        
        u=coef*randn(n,v)*chrho;
        covu=cov(u);
        
        [OMG,~] = eig(covu);
        OMG3D(:,:,j)=OMG;
    end
else
    for j=1:k
        gam=sort(rand(v,1));
        gam=gam./(prod(gam)^(1/v));
        GAM3D(:,:,j)=diag(gam);
        % GAM2D(:,j)=gam;
        
        u=coef*randn(n,v)*chrho;
        covu=cov(u);
        
        [OMG,~] = eig(covu);
        OMG3D(:,:,j)=OMG;
    end
end

% Extract modeltype
modeltype=pa.pars;

if strcmp(modeltype,'VVE')
    
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Equal orientation
        OMG3Dj=OMG3D(:,:,1);
        % S(:,:,j)=(lmdj^(1/v)) *  OMG3Dj * GAM3Dj * OMG3Dj'; TODOTODO
        S(:,:,j)=lmdj *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'EVE')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Equal orientation
        OMG3Dj=OMG3D(:,:,1);
        S(:,:,j)=lmdj *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'VVV')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Variable orientation
        OMG3Dj=OMG3D(:,:,j);
        S(:,:,j)=lmdj *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'EVV')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Variable orientation
        OMG3Dj=OMG3D(:,:,j);
        S(:,:,j)=lmdj *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'VEE')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Equal orientation
        OMG3Dj=OMG3D(:,:,1);
        S(:,:,j)=lmdj *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'EEE')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Equal orientation
        OMG3Dj=OMG3D(:,:,1);
        S(:,:,j)=lmdj *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'VEV')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Variable orientation
        OMG3Dj=OMG3D(:,:,j);
        S(:,:,j)=lmdj *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'EEV')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Variable orientation
        OMG3Dj=OMG3D(:,:,j);
        S(:,:,j)=lmdj *  OMG3Dj * GAM3Dj * OMG3Dj';
    end
    
elseif    strcmp(modeltype,'VVI')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Axis-Aligned
        S(:,:,j)=lmdj *   GAM3Dj ;
    end
    
elseif    strcmp(modeltype,'EVI')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Variable shape
        GAM3Dj=GAM3D(:,:,j);
        % Axis-Aligned
        S(:,:,j)=lmdj *   GAM3Dj ;
    end
    
elseif    strcmp(modeltype,'VEI')
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Axis-Aligned
        S(:,:,j)=lmdj *   GAM3Dj ;
    end
    
elseif    strcmp(modeltype,'EEI')
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Equal shape
        GAM3Dj=GAM3D(:,:,1);
        % Axis-Aligned
        S(:,:,j)=lmdj *   GAM3Dj ;
    end
    
elseif    strcmp(modeltype,'VII')
    eyev=eye(v);
    for j=1:k
        % Variable determinant (volume)
        lmdj=lmd(j);
        % Spherical shape
        S(:,:,j)=lmdj *   eyev ;
    end
    
elseif    strcmp(modeltype,'EII')
    eyev=eye(v);
    for j=1:k
        % Equal determinant (volume)
        lmdj=lmd(1);
        % Spherical shape
        S(:,:,j)=lmdj *   eyev ;
    end
else
    error('FSDA:genSigmaGPCM:WrongModelSupplied','Model supplied id not recognized');
end

end
%FScategory:CLUS-RobClaMULT