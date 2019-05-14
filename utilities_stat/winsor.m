function Xwins = winsor(X, p, dim)
%winsor returns a winsorized copy of input 
%
% Measurement data often contains "outliers," sample points rather far
% outside the range containing the majority of the points. While expected
% both from theory and experience, these outliers, for small or
% medium-sized samples, tend to distort statistical data such as the mean
% value. One of the standard methods dealing with this problem for (real)
% continuous scales is clamping the outliers. Function winsor sets all
% data points below or above a given quantile to these quantiles. (This
% operation is named after its inventor, Charles P. Winsor.)
%
%<a href="matlab: docsearchFS('winsor')">Link to the help function</a>
%
% Required input arguments:
%
% 
%    X  :   Input array. Vector | matrix | 3D array.
%           Input array, specified as a vector, matrix, or 3D  array.
%           Data Types -single | double | int8 | int16 | int32 | int64 |uint8 | uint16 | uint32 | uint64
%    p  :   Percentiles to use to winsorize. Vector of length 2 with elements in the interval [0 100].
%           Vector of length(2) containing cut-off percentiles (left,
%           right). Note that p(2) must be greater than p(1).
%           Data Types - single | double | int8 | int16 | int32 | int64 |uint8 | uint16 | uint32 | uint64
% 
% Optional input arguments:
%
%  dim  :   Dimension to operate along. Positive integer scalar. 
%           If no value is specified, then the default is the first
%           array dimension whose size does not equal 1. In other words, if
%           the input is a row (column) vector it winsorizes this row
%           (column) vector. If the input is a nxp matrix it winsorizes the
%           columns. If the input is a 3D arrary of size nxpxq it
%           winsorizes each column of each slice of the 3D array.
%           Example - 2 
%           Data Types -single | double | int8 | int16 | int32 | int64 |uint8 | uint16 | uint32 | uint64
%
% Output:
%
% Xwins  :  X winsorized. Array of the same size of X. In array Xwins all
%           entries smaller than the p(1) quantile have been replaced by
%           this value and likewise for all entries larger than the p(2)
%           quantile. Row, or columns are winsorized depending on optional
%           input scalar dim.
%
%
% See also: trimmean
%
% References:
%
% Maronna, R.A., Martin D. and Yohai V.J. (2006), "Robust Statistics, Theory
% and Methods", Wiley, New York.
%
% Copyright 2008-2019.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('winsor')">Link to the help function</a>
%
%$LastChangedDate::                      $: Date of the last commit
%

% Examples:

%{ 
    %% Winsorize each column of matrix X.
    rng('default') % Reinitialize the random number generator to its startup configuration
    rng(100)
    X = randn(100,5);
    % Contaminate 5th and 10th row of matrix X
    X(5,:)=5;
    X(10,:)=-6;
    p=[1 99];
    Y = winsor(X,p);
    subplot(2,1,1)
    boxplot(X)
    title('Original data')
    subplot(2,1,2)
    boxplot(Y)
    title('Winsorized data')
%}

%{ 
    % Winsorize each row of 3D array X.
    rng(100)
    X = randn(20,100,3);
    Y = winsor(X,p,2);
    disp('Number of unique elements in second rows of first slice of 3D array X')
    disp(length(unique(Y(2,:,1)')))
%}

%% Beginning of code

if nargin < 2
   error('FSDA:winsor:InvalidArg','Input argument "p" is undefined')
end 
if ~isvector(p)
   error('FSDA:winsor:InvalidArg','Input argument "p" must be a vector')
end  
if p(1) < 0 || p(1) > 100
   error('FSDA:winsor:InvalidArg','Left cut-off percentile is out of [0,100] range')
end  
if p(2) < 0 || p(2) > 100
   error('FSDA:winsor:InvalidArg','Right cut-off percentile is out of [0,100] range')
end  
if p(1) > p(2)
   error('FSDA:winsor:InvalidArg','Left cut-off percentile exceeds right cut-off percentile')
end

if nargin == 2 && isvector(X) % Input is a  vector
    Xwins=winsorcore(X,p);
else % Input is at least a two dimensional array
    if nargin == 2    % Determine first nonsingleton dimension
        dim = find(size(X)~=1,1);
    end
    s = size(X);
    if dim > length(s)  % If dimension is too high, just return input.
        Xwins = X;
        return
    end
    if length(s)==2
        if dim==1  % Input is a matrix  dim=1 compute along columns
            Xwins=zeros(s(1),s(2));
            for j=1:s(2)
                Xwins(:,j)=winsorcore(X(:,j),p);
            end
        else  % Input is a matrix  dim=2 compute along rows
            Xwins=zeros(s(1),1);
            for i=1:s(1)
                Xwins(i,:)=winsorcore(X(i,:),p);
            end
        end
    elseif length(s)==3
            Xwins=zeros(s);

        if dim==1 % Input is a 3D array dim=1 compute along columns
            for k=1:s(3)
                for j=1:s(2)
                    Xwins(:,j,k)=winsorcore(X(:,j,k),p);
                end
            end
            
        elseif dim==2  % Input is a 3D array dim=2 compute along rows
            for k=1:s(3)
                for i=1:s(1)
                    Xwins(i,:,k)=winsorcore(X(i,:,k),p);
                end
            end
            
        else % Input is a 3D array dim=3 compute along 3rd dim
            
            for i=1:s(1)
                for j=1:s(2)
                    Xwins(i,j,:)=winsorcore(X(i,j,:),p);
                end
            end
        end
    else
        error('FSDA:winsor:WrongInput','Not implemented for array of size greater than 3')
    end
end
end

% winsorcore is the function which effectively performs the winsorization
% of vector X
function xwins=winsorcore(X,p)
p = prctile(X,p);
i1 = X < p(1); v1 = min(X(~i1));
i2 = X > p(2); v2 = max(X(~i2));
xwins = X;
xwins(i1) = v1;
xwins(i2) = v2;

end

%FScategory:UTISTAT