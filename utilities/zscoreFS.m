function [Z,mu,sigma] = zscoreFS(X,loc,scale,dim)
%zscoreFS computes (robust) standardized z scores
%
%<a href="matlab: docsearchFS('zscorefs')">Link to the help function</a>
%
% Required input arguments:
%
%    X: vector of length(n) or data matrix containing n observations on v
%    variables or 3D array of size n-by-v-by-r.
%       Missing values (NaN's) and infinite values (Inf's) are allowed,
%       since observations (rows) with missing or infinite values will
%       automatically be excluded from the computations.
%
%   Z = zscoreFS(X) returns a centered, scaled version of X, with the same size
%   as X. For vector input X, Z is the vector of z-scores
%
%      (X-median(X)) ./ mad(X).
%
%   Z=zscoreFS(X,loc,scale) returns a centered, scaled version of X, the
%   same size as X using location and scale are specified in input
%   parameters 'loc' and 'scale'. For vector input X, Z is the vector of
%   z-scores
%
%      (X-location(X)) ./ scale(X).
%
%    loc is a string which specifies the location measure to use
%     'median' (default).
%     'mean'
%    scale is a string which specifies the dispersion measure to use
%     'mad' (default). Traditional mad is Me(|x_i-Me(X)|)/norminv(3/4)
%     'Qn' first quartile of interpoint distances |x_i-x_j|. See function Qn.m  
%     'Sn' robust Gini's average difference index. See function Sn.m
%     'std' Unbiased standard deviations. See function std.m ()
%     'modmadp'. Modified mad where the last letter(s) p of string modmap
%                 is (are) a number converted to string necessary to
%                 compute the modified MAD. 
%       Modified MAD = (order statistic ceil((n+p-1)/2) of |x_i-Me(X)|
%                 + order statistic floor((n+p-1)/2+1) of |x_i-Me(X)|)
%                 / (2* ?) where ?= norminv(0.5*((n+p-1)/(2*n)+1))
%                  Note that p is different from v (columns of X if X is a
%                  matrix) and must be supplied by the user.
%       For example if p=5 then the user can supply the string 'modmad5'
%        as follows.  
%             p=5; modmadp=['modmap' num2str(p)];
%
%
%   Z=zscoreFS(X,loc,scale) computes robust standardized zscores using the
%   estimates of location and scale specified in loc and scale strings. If
%   X is a 2D matrix, zscores are computed using loc and scale along each
%   column of X. If X is a 3D array zscores are
%   computed using the location and scale along the first
%   non-singleton dimension. For example if X is n-by-v-by-r (with n>1) and
%   loc='median'; n-by-r medians are computed for each of the n rows of X
%   and each third dimension r
%
%
%   Z=zscoreFS(X,loc) computes standardized zscores using the
%   estimates of location specified in loc and the mad as measure of
%   dispersion
%
%   Z=zscoreFS(X,[],scale) computes standardized zscores using as
%   estimates of location the medians and the measure of dispersion
%   specified in scale
%
%   [Z,mu,sigma] = zscoreFS(X) also returns median(X) in mu and mad in
%   sigma
%
%   [Z,mu,sigma] = zscoreFS(X,loc,scale) also returns the estimates of location
%   in mu and of scale in sigma as specified in loc and scale strings
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
%   scale
%
%
% See also zscore, MCD, Smult, MMmult, FSM
%
%
% Copyright 2008-2015.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('zscorefs')">Link to the help function</a>
% Last modified 06-Feb-2015
%

% Examples

%{
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
    % zscoreFS (remove the medians and
    % divide by Qn)
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
    % Examples with 3D arrays
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
    % 3D arrays with dim=1, dim=2 and dim=3
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
    % Example of use of modmad as a scale measure
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
            x = sort(x);
            locest = x(half+1,:);
            if 2*half == n       % Average if even number of elements
                locest =(x(half,:)+locest)/2;
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
                x = sort(x);
                locestmean = x(half+1,:);
                if 2*half == n       % Average if even number of elements
                    locestmean =(x(half,:)+locestmean)/2;
                end
                signres = bsxfun(@minus,x, locestmean);
            else
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