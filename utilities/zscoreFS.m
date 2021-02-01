function [Z,mu,sigma] = zscoreFS(X,loc,scale,dim)
%zscoreFS computes (robust) standardized z scores
%
%<a href="matlab: docsearchFS('zscoreFS')">Link to the help function</a>
%
%
%    X can be a vector of length(n) or data matrix containing n observations on v
%       variables or 3D array of size n-by-v-by-r.
%   Z = zscoreFS(X) returns a centered, scaled version of X, with the same size
%   as X. For vector input X, Z is the vector of z-scores
%
%      (X-median(X)) ./ (1.4826* mad(X)).
%
%   Z=zscoreFS(X,loc,scale) returns a centered, scaled version of X, the
%   same size as X using location and scale are specified in input
%   parameters 'loc' and 'scale'. For vector input X, Z is the vector of
%   z-scores
%
%      (X-location(X)) ./ scale(X).
%
%   where scaled(X) is the corrected estimator of scale (corrected in the
%   sense that it is multiplied by a coefficient to achieve consistency for 
%   normally distributed data).  
%
%
%   Z=zscoreFS(X,loc,scale) computes robust standardized zscores using the
%   estimates of location and scale specified in loc and scale strings. If
%   X is a 2D matrix, zscores are computed using loc and scale along each
%   column of X. If X is a 3D array zscores are
%   computed using the location and scale along the first
%   non-singleton dimension. For example if X is n-by-v-by-r (with n>1) and
%   loc='median'; n-by-r medians are computed for each of the n rows of X
%   and each third dimension r.
%
%
%   Z=zscoreFS(X,loc) computes standardized zscores using the
%   estimates of location specified in loc and the mad as measure of
%   dispersion.
%
%
%   [Z,mu,sigma] = zscoreFS(X) also returns median(X) in mu and mad in
%   sigma.
%
%   [Z,mu,sigma] = zscoreFS(X,loc,scale) also returns the estimates of location
%   in mu and of scale in sigma as specified in loc and scale strings.
%
%   Z=zscoreFS(X,loc,scale,dim) computes robust standardized zscores along
%   the dimension dim of X using the estimates of location and scale
%   specified in loc and scale strings. dim standardizes X by working along
%   the dimension dim of X. For example if X is a two dimensional matrix
%   dim=2 (default) standardizes the columns of X else if dim=1
%   standardizes the rows. If X is a three dimensional dim = 1 standardizes
%   the columns, dim =2 standardizes the rows and dim =3 standardizes the
%   third dimension.
%
%   zscoreFS is an extension of function zscore of statistic toolbox
%   because it enables to specify alternative measures of location and
%   scale.
%
%
%  Required input arguments:
%  
% X :           Input data. Vector or Matrix or 3D array. Vector  of
%               length n or data matrix containing n
%               observations on v variables or 3D array of size
%               n-by-v-by-r.
%               Missing values (NaN's) and infinite values (Inf's) are
%               allowed, since observations (rows) with missing or infinite
%               values will automatically be excluded from the
%               computations.
%                Data Types - single|double
%
%
%  Optional input arguments:
%
%   loc : location measure to use. 'median' (default) or 'mean'.
%         String which specifies the location measure to use. The default
%         value is 'median'. 
%               Example - 'median'
%               Data Types - character
% scale : scale measure to use. 'mad' (default) or 'Qn' or 'Sn' or 'std' or
%         moddmadp'.
%         String which specifies the dispersion measure to use
%           'mad' is the default. Traditional (corrected) mad is
%           $Me(|x_i-Me(X)|)/norminv(3/4)$;
%           'Qn' first quartile of interpoint distances $|x_i-x_j|$ corrected
%           by the consistency factor. See function Qn.m;
%           'Sn' robust Gini's average difference index corrected by the
%           consistency factor. See function Sn.m;
%           'std' Unbiased standard deviations. See function std.m; 
%           'modmadp'. Modified mad where the last letter(s) p of string modmap
%                 is (are) a number converted to string necessary to
%                 compute the modified MAD. 
%       Modified MAD = (order statistic $ceil((n+p-1)/2)$ of $|x_i-Me(X)|$
%                 + order statistic $floor((n+p-1)/2+1)$ of $|x_i-Me(X)|)$
%                 / $(2 \sigma)$ where $\sigma=
%                 norminv(0.5*((n+p-1)/(2*n)+1))$.
%                  Note that $p$ is different from $v$ (columns of X if X is a
%                  matrix) and must be supplied by the user.
%                   For example if p=5 then the user can supply the string 'modmad5'
%                   as follows.  p=5; modmadp=['modmap' num2str(p)];
%               Example - 'mad'
%               Data Types - character
%  dim  :   Dimension to operate along. Positive integer scalar. 
%           Dimension to operate along, specified as a positive integer
%           scalar. If no value is specified, then the default is the first
%           array dimension whose size does not equal 1.
%               Example - 2
%           Data Types -single | double | int8 | int16 | int32 | int64 |uint8 | uint16 | uint32 | uint64
%
%
%  Output: 
%
%       Z : centered, scaled version of X. Array with the same dimension as input X.
%           Array with the same size as X using location and scale are specified in input
%           parameters 'loc' and 'scale'. For vector input X, Z is the vector of
%            z-scores
%           (X-location(X)) ./ scale(X).
%   mu : location estimate. Scalar, vector or matrix depending on the size of input matrix X.
%           Estimates of location specified in loc input string.
%  sigma : scale estimate. Scalar, vector or matrix depending on the size of input matrix X.
%           Estimates of scale specified in scale input string.
%
% See also: zscore, Qn, Sn, MCD, Smult, MMmult, FSM
%
%
% References:
%
%
%
% Copyright 2008-2021.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('zscoreFS')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples



%{
    % Scale using medians and mads.
    % zscoreFS with all default options (that is remove the medians and
    % divide by mads)
    n=200;
    v=3;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+10;
    [out]=zscoreFS(Ycont);
%}

%{
    % Scale using mean and mads.
    % Computes standardized zscores using mean and mads
    % estimates of location the medians and the measure of dispersion
    % specified in scale
    loc='mean'
    X=randn(10,2);
    Z=zscoreFS(X,loc,'mad'); 
%}

%{
    % Remove the medians and divide by Qn.
    n=200;
    v=1;
    randn('state', 123456);
    Y=randn(n,v);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:)=Ycont(1:5,:)+10;
    [out]=zscoreFS(Ycont,[],'Qn');
    % Alternatively it is possible to use the following sintax
    [out]=zscoreFS(Ycont,'median','Qn');
%}

%{
    % Examples with 3D arrays.
    n=200;
    v=3;
    q=5;
    randn('state', 123456);
    Y=randn(n,v,q);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:,:)=Ycont(1:5,:,:)+10;
    [out1,Mu,Sigma]=zscoreFS(Ycont,[],'Sn',1);
    % [out,Mu1,Sigma1]=zscoreFS(Ycont,[],'Sn',1);
%}

%{
    % Report also location and scale measures which have have been used.
    % zscoreFS produces the same output as function zscore of statistics
    % toolbox if centroid is arithmetic mean and scale measure is the
    % standard deviation
    X=randn(10,3,6);
    [Z,mu,sig]=zscoreFS(X,'mean','std',3);
    [Z1,mu1,sig1]=zscore(X,[],3);
    if isequal(Z,Z1) + isequal(mu,mu1) + isequal(sig,sig) ==3
        disp('Everything is equal')
    else
        disp('Equality not reached')
    end
%}

%{
    % 3D arrays with dim=1, dim=2 and dim=3.
    n=200;
    v=3;
    q=5;
    randn('state', 123456);
    Y=randn(n,v,q);
    % Contaminated data
    Ycont=Y;
    Ycont(1:5,:,:)=Ycont(1:5,:,:)+10;
    scale='Qn';
    loc='mean';
    dim=2; % work along rows
    [Z,Mu1,Sigma1]=zscoreFS(Ycont,loc,scale,dim);
    isequal(Z(3,:,2)',zscoreFS(Ycont(3,:,2),loc,scale))

    scale='Qn';
    loc='median';
    dim=1; % work along columns
    [Z,Mu1,Sigma1]=zscoreFS(Ycont,loc,scale,dim);
    isequal(Z(:,2,4),zscoreFS(Ycont(:,2,4),loc,scale))

    scale='Sn';
    loc='median';
    dim=3; % work along third dimension
    [Z,Mu1,Sigma1]=zscoreFS(Ycont,loc,scale,dim);
    isequal(squeeze(Z(7,2,:)),zscoreFS(squeeze(Ycont(7,2,:)),loc,scale))
%}

%{
    % Example of use of modmad as a scale measure.
    p=3;
    X=randn(100,p);
    loc='median';
    scale=['modmad' num2str(p)];
    % Project the data using v vectors   
    v=10;
    proj=randn(p,v);
    Y=X*proj;
    % Standardize the n projected points using median and modified MAD
    % Note that Y has v columns but the original matrix X has p columns
    [Z,Mu1,Sigma1]=zscoreFS(Y,loc,scale);   

%}


%% Beginning of code

if nargin <2
    % location and scale not specified. Use median and mad
    loc='median';
    scale='mad';
elseif nargin ==2
    scale='mad';
elseif nargin>=3
    if isempty(loc)
        loc='median';
    end
end

if nargin < 4 && iscolumn(X) % Input is a column  vector
    [Z,mu,sigma]=zscoreFScore(X,loc,scale);
elseif nargin < 4 && isrow(X) % Input is a row vector
    [Z,mu,sigma]=zscoreFScore(X',loc,scale);
else % Input is at least a two dimensional array
    if nargin < 4    % Determine first nonsingleton dimension
        dim = find(size(X)~=1,1);
    end
    s = size(X);
    if dim > length(s)  % If dimension is too high, just return input.
        Z = X;
        return
    end
    if length(s)==2
        if dim==1  % Input is a matrix  dim=1 compute along columns
            [Z,mu,sigma]=zscoreFScore(X,loc,scale);
        else  % Input is a matrix  dim=2 compute along rows
            [Z,mu,sigma]=zscoreFScore(X',loc,scale);
        end
    elseif length(s)==3
        if dim==1 % Input is a 3D array dim=1 compute along columns
            Z=zeros(s);
            mu=zeros(1,s(2),s(3)); % first dimension is collapsed
            sigma=mu;
            for k=1:s(3)
                [Zk,muk,sigmak]=zscoreFScore(X(:,:,k),loc,scale);
                Z(:,:,k)=Zk;
                mu(:,:,k)=muk;
                sigma(:,:,k)=sigmak;
            end
            
        elseif dim==2  % Input is a 3D array dim=2 computes along rows
            Z=zeros(s);
            mu=zeros(s(1),1,s(3)); % second dimension is collapsed
            sigma=mu;
            
            for k=1:s(3)
                [Zk,muk,sigmak]=zscoreFScore(X(:,:,k)',loc,scale);
                Z(:,:,k)=Zk';
                mu(:,:,k)=muk';
                sigma(:,:,k)=sigmak';
            end
            
        else % Input is a 3D array dim=3 compute along 3rd dim
            Z=zeros(s);
            mu=zeros(s(1),s(2),1); % third dimension is collapsed
            sigma=mu;
            
            for i=1:s(1)
                [Zi,mui,sigmai]=zscoreFScore(squeeze(X(i,:,:))',loc,scale);
                Z(i,:,:)=Zi';
                mu(i,:)=mui;
                sigma(i,:)=sigmai;
            end
        end
    else
        error('FSDA:zscoreFS:WrongInput','Not implemented for array of size greater than 3')
    end
end

    function [z,locest,scaleest]=zscoreFScore(x,loc,scale)
        n=size(x,1);
        half = floor(n/2);
        
        if strcmp(loc,'median')
            xsor = sort(x);
            locest = xsor(half+1,:);
            if 2*half == n       % Average if even number of elements
                locest =(xsor(half,:)+locest)/2;
            end
        elseif strcmp(loc,'mean')
            locest =sum(x,1)/n;
        end
        
        madscale=strcmp(scale,'mad');
        if length(scale)>3
            modmadscale=strcmp(scale(1:6),'modmad');
        else
            modmadscale=0;
        end
        if  madscale || modmadscale
            fln2=floor(n/2);
            if modmadscale
                p=str2double(scale(7:end));
                if fln2>p
                    n1 = ceil((n+p-1)/2);
                    n2 = floor((n+p-1)/2)+1;
                    nb = (n+p-1)/2;
                    % beta = constant necessary to rescale the modified MAD
                    beta = norminv(0.5*(nb/n+1),0,1);
                else
                    disp('Modified MAD has been selected of a matrix which has a number of columns greater than half the number of rows')
                    disp('Use simple MAD instead')
                    error('FSDA:zscoreFS:WrongInput','Stop execution')
                end
                
            else % Simple MAD has been selected
                nx = fln2;
                if 2*nx == n       % Average if even number of elements
                    n1 =nx;
                    n2=nx+1;
                else
                    n1=nx+1;
                    n2=nx+1;
                end
                beta = norminv(3/4,0,1);
            end
            
            if strcmp(loc,'mean')
              % In this case the median has not been calculated yet
                xsor = sort(x);
                locestmedian = xsor(half+1,:);
                if 2*half == n       % Average if even number of elements
                    locestmedian =(xsor(half,:)+locestmedian)/2;
                end
                signres = bsxfun(@minus,x, locestmedian);
            else
                % in this case the median had already been calculated and
                % put in variable locest.
                signres = bsxfun(@minus,x, locest);
            end
            
            % MADs or modified MADs
            res = abs(signres);
            ress = sort(res);
            % Divide by beta (asymptotic consistency factor for (modified) MAD)
            scaleest = (ress(n1,:)+ress(n2,:))/(2*beta);
            
        elseif strcmp(scale,'Sn')
            scaleest=Sn(x,1);
            
        elseif strcmp(scale,'Qn')
            scaleest=Qn(x,1);
            
        elseif strcmp(scale,'std')
            scaleest=std(x);
        else
            error('FSDA:zscoreFS:WrongScale','Scale estimate can be Sn, Qn or std')
        end
        
        z = bsxfun(@minus,x, locest);
        z = bsxfun(@rdivide, z, scaleest);
    end
end
%FScategory:UTIGEN