function y=Qn(X,dim)
%Qn robust estimator of scale (first quartile of interpoint distances $|x_i-x_j|$)
%
%<a href="matlab: docsearchFS('Qn')">Link to the help function</a>
%
%
% Required input arguments:
% 
%    X  :   Input array. Vector | matrix | 3D array.
%           Input array, specified as a vector, matrix, or 3D  array.
%           For vectors, Qn(X) is the scale estimator of the elements in X. For
%           matrices, Qn(X) is a row vector containing the scale estimator
%           value of each column.  For 3D arrays, Qn(X) is the robust scale
%           estimator of the elements along the first non-singleton dimension
%           of X.
%           Data Types -single | double | int8 | int16 | int32 | int64 |uint8 | uint16 | uint32 | uint64
%  dim  :   Dimension to operate along. Positive integer scalar. 
%           Dimension to operate along, specified as a positive integer
%           scalar. If no value is specified, then the default is the first
%           array dimension whose size does not equal 1.
%           Data Types -single | double | int8 | int16 | int32 | int64 |uint8 | uint16 | uint32 | uint64
% 
% Optional input arguments:
%
% Output:
%
%    y  :    robust estimator of scale. Scalar | Vector or 3D array. 
%            Qn(X,dim) takes the robust estimator of scale along the
%            dimension dim of X.
%
% More About:
%
%   $Q_n$ is the first quartile of the distances { $|x_i-x_j|$; $i <j$} Note that
%   $Q_n$ does not need any location estimate. More in detail, let 
%   $d_{(1)} \leq d_{(2)} \leq ... \leq d_{(m)}$ the ordered values of the
%   $m$
%   differences $|x_i-x_j|$ with $i>j$ and $m = {n \choose 2}$. $Q_n=d_{(k)}$ where
%   $k= {[n/2]+1 \choose 2}$. Since $k$ is approximately $m/4$, $Q_n$ is approximately
%   the first quartile of the ordered distances $d_{(1)} \leq d_{(2)} \leq
%   ... \leq d_{(m)}$. $Q_n$ is multiplyed by $c$ and $c_n$. 
%   $c$ is the so called
%   asymptotic consistency factor and is equal to 2.2219 while $c_n$ is a
%   finite sample correction factor to make the estimator unbiased.
%
%
% See also: Sn
%
% References:
%
% Rousseeuw P.J. and Croux C., (1993), Alternatives to the median absolute deviation,
% Journal of American Statistical Association 88, 1273-1283
% Croux C. and Rousseeuw P.J.(1992)  Time-efficient algorithms for two
% highly robust estimators of scale, in Computational Statistics, Volume 1,
% eds. Y . Dodge and J. Whittaker, Heidelberg: Physika-Verlag, 41 1-428.
%
% Copyright 2008-2016.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('Qn')">Link to the help function</a>
% Last modified mer 25 mag 2016 18:19:58
%

% Examples:

%{ 
    %% Qn with all default otpions.
    X = [1 2 4 4  7;
         3 4 6 6  8;
         5 6 8 8  10;
         5 7 10 12 1500];
    y1=Qn(X)
    y2=Qn(X,2)
%}


%% Beginning of code

if nargin == 1 && iscolumn(X) % Input is a column  vector
    y=Qncore(X);
elseif nargin == 1 && isrow(X) % Input is a row vector
    y=Qncore(X');
else % Input is at least a two dimensional array
    if nargin == 1    % Determine first nonsingleton dimension
        dim = find(size(X)~=1,1);
    end
    s = size(X);
    if dim > length(s)  % If dimension is too high, just return input.
        y = X;
        return
    end
    if length(s)==2
        if dim==1  % Input is a matrix  dim=1 compute along columns
            y=zeros(1,s(2));
            for j=1:s(2)
                y(:,j)=Qncore(X(:,j));
            end
        else  % Input is a matrix  dim=2 compute along rows
            y=zeros(s(1),1);
            for i=1:s(1)
                y(i)=Qncore(X(i,:)');
            end
        end
    elseif length(s)==3
        if dim==1 % Input is a 3D array dim=1 compute along columns
            y=zeros(1,s(2),s(3));
            for k=1:s(3)
                for j=1:s(2)
                    y(1,j,k)=Qncore(X(:,j,k));
                end
            end
            
        elseif dim==2  % Input is a 3D array dim=2 compute along rows
            y=zeros(s(1),1,s(3));
            for k=1:s(3)
                for i=1:s(1)
                    y(i,1,k)=Qncore(X(i,:,k)');
                end
            end
            
        else % Input is a 3D array dim=3 compute along 3rd dim
            y=zeros(s(1),s(2));
            
            for i=1:s(1)
                for j=1:s(2)
                    y(i,j)=Qncore(X(i,j,:)');
                end
            end
        end
    else
        error('FSDA:Qn:WrongInput','Not implemented for array of size greater than 3')
    end
end


    function s=Qncore(x)
        
        n=length(x);
        % Do binning for big n
        if n>10000
            sy=sort(x);
            nbins=floor(n/10);
            xbinned=zeros(nbins,1);
            ninbins=floor(n/nbins);
            for ii=1:nbins
                if (mod(n,nbins)~=0 && ii==nbins)
                    xbinned(ii)=median(sy((ii-1)*ninbins+1:n));
                else
                    xbinned(ii)=median(sy((ii-1)*ninbins+1:ii*ninbins));
                end
            end
            x=xbinned; % Redefine x with binned x
            n=nbins; % Redefine n with number of bins
        end
        h=floor(n/2)+1;
        kk=0.5*h*(h-1);
        
        % Compute the n*(n-1)/2 pairwise ordered distances
        % Use function pdist of statistics toolbox
        distord=sort(pdist(x,'cityblock'));
        
        %        If statistic toolbox is not present it is possible to use the following code
        %         distord = zeros(1,n*(n-1)./2);
        %         kkk = 1;
        %         for iii = 1:n-1
        %             d = abs(x(iii) - x((iii+1):n));
        %             distord(kkk:(kkk+n-iii-1)) = d;
        %             kkk = kkk + (n-iii);
        %         end
        %         distord=sort(distord);
        
        s=2.2219*(distord(kk));

        % Multiply the estimator also by cn a finite sample correction
        % factor to make the estimator unbiased for finite samples (see p. 10
        % of Croux and Rousseeuw, 1992) or
        % http://www.google.it/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CDUQFjAA&url=http%3A%2F%2Fwww.researchgate.net%2Fpublication%2F228595593_Time-efficient_algorithms_for_two_highly_robust_estimators_of_scale%2Ffile%2F79e4150f52c2fcabb0.pdf&ei=ZCE5U_qHIqjU4QTMuIHwAQ&usg=AFQjCNERh4HiLgtkUGF1w4JU1380xhvKhA&bvm=bv.63808443,d.bGE

        switch n
            case 1
                error('FSDA:Qn:TooSmalln','Sample size too small');
            case 2
                dn=0.399;
            case 3
                dn=0.994;
            case 4
                dn=0.512;
            case 5
                dn=0.844;
            case 6
                dn=0.611;
            case 7
                dn=0.857;
            case 8
                dn=0.669;
            case 9
                dn=0.872;
                
            otherwise
                if (mod(n,2)==1)
                    dn=n/(n+1.4);
                elseif (mod(n,2)==0)
                    dn=n/(n+3.8);
                end
        end
        s=dn*s;
    end
end
%FScategory:UTISTAT