function y=Sn(X,dim)
%Sn robust estimator of scale (robust version of Gini's average difference)
%
%<a href="matlab: docsearchFS('Sn')">Link to the help function</a>
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
%            Sn(X,dim) takes the robust estimator of scale along the
%            dimension dim of X.
%
% More About:
%
%   For vectors, Sn(X) is the scale estimator of the elements in X. For
%   matrices, Sn(X) is a row vector containing the scale estimator value of
%   each column.  For 3D arrays, Sn(X) is the robust scale estimator of the
%   elements along the first non-singleton dimension of X.
%
%   Sn(X,dim) takes the robust estimator along the dimension dim of X.
%
%   $Sn= cn \times c \times med_i { med_j |x_i-x_j|}$, $i=1,2, ...n$, $j=1, 2, ..., n$.
%   For each $i$ we compute the median of $|x_i-x_j|$, $j=1, 2, ..., n$.
%   This yields $n$ numbers, the median of which gives the final estimate of
%   $S_n$. This estimator is the robust version of Gini's average difference,
%   which one would obtain when replacing medians by averages
%   More in detail $Sn = cn \times c \times lomed_i { highmed_j
%   |x_i-x_j|}$, $i=1,2, ...n$, $j=1, 2, ..., n$, where $lomed$ (low
%   median) is $[(n+1)/2]$-th order statistic out of $n$ numbers) and
%   $himed$ (high median) is the $([n/2]+1)$-th order statistic out of the
%   $n$ numbers. $c$ is the so called asymptotic consistency factor and is
%   equal to 1.1926 while $cn$ is a finite sample correction factor to make
%   the estimator unbiased.
% 
% See also: Qn
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
%<a href="matlab: docsearchFS('Sn')">Link to the help function</a>
% Last modified mer 25 mag 2016 18:19:58
%

% Examples:

%{ 
    %% Sn with all default otpions.
    X = [1 2 4 4  7;
         3 4 6 6  8;
         5 6 8 8  10;
         5 7 10 12 1500];
    y1=Sn(X)
    y2=Sn(X,2)
%}


%% Beginning of code

if nargin == 1 && iscolumn(X) % Input is a column  vector
    y=Sncore(X);
elseif nargin == 1 && isrow(X) % Input is a row vector
    y=Sncore(X');
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
                y(:,j)=Sncore(X(:,j));
            end
        else  % Input is a matrix  dim=2 compute along rows
            y=zeros(s(1),1);
            for i=1:s(1)
                y(i)=Sncore(X(i,:)');
            end
        end
    elseif length(s)==3
        if dim==1 % Input is a 3D array dim=1 compute along columns
            y=zeros(1,s(2),s(3));
            for k=1:s(3)
                for j=1:s(2)
                    y(1,j,k)=Sncore(X(:,j,k));
                end
            end
            
        elseif dim==2  % Input is a 3D array dim=2 compute along rows
            y=zeros(s(1),1,s(3));
            for k=1:s(3)
                for i=1:s(1)
                    y(i,1,k)=Sncore(X(i,:,k)');
                end
            end
            
        else % Input is a 3D array dim=3 compute along 3rd dim
            y=zeros(s(1),s(2));
            
            for i=1:s(1)
                for j=1:s(2)
                    y(i,j)=Sncore(X(i,j,:)');
                end
            end
        end
    else
        error('FSDA:Sn:WrongInput','Not implemented for array of size greater than 3')
    end
end


    function s=Sncore(x)
        
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
        
        dist_xi_xj_sor=sort(abs(bsxfun(@minus,x,x')));
        
        % For each i compute the median of |x_i-x_j| that is take di order
        % statistic of rank [(n+1)/2]
        med_j=dist_xi_xj_sor(floor((n+1)/2),:);
        
        med_i=sort(med_j);
        s=1.1926*(med_i(floor(n/2)+1));
        
        % Multiply the estimator also by cn a finite sample correction
        % factor to make the estimator unbiased for finite samples (see p. 3
        % of Croux and Rousseeuw, 1992) or
        % http://www.google.it/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0CDUQFjAA&url=http%3A%2F%2Fwww.researchgate.net%2Fpublication%2F228595593_Time-efficient_algorithms_for_two_highly_robust_estimators_of_scale%2Ffile%2F79e4150f52c2fcabb0.pdf&ei=ZCE5U_qHIqjU4QTMuIHwAQ&usg=AFQjCNERh4HiLgtkUGF1w4JU1380xhvKhA&bvm=bv.63808443,d.bGE

        cn=1;
        switch n
            case 2
                cn=0.743;
            case 3
                cn=1.851;
            case 4
                cn=0.954;
            case 5
                cn=1.351;
            case 6
                cn=0.993;
            case 7
                cn=1.198;
            case 8
                cn=1.005;
            case 9
                cn=1.131;
            otherwise
                if (mod(n,2)==1)
                    cn=n/(n-0.9);
                end
        end
        s=cn*s;
    end
end
%FScategory:UTISTAT