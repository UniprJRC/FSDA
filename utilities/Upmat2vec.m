function y=Upmat2vec(A)
%Upmat2vec extracts in a vector, linear indexes or elements above main diagonal of a square matrix
%
%
%<a href="matlab: docsearchFS('upmat2vec')">Link to the help function</a>
%
%  Required input arguments:
%
%         A   : scalar, or square matrix  of order k (k>1).
%
%               If A is a scalar, A is equal to the number of rows
%               (columns) of a matrix whose linear indexes of the elements
%               above diagonal have to be found. In the case the procedure
%               returns the (A*(A-1)/2) linear indexes of the elements
%               above diagonal.
%
%               If A is a square matrix of order k (k>1), the procedure
%               returns the elements above diagonal of matrix A in a vector
%               of length (k*(k-1)/2);
%
%  Output:
%
%        y  : vector containing the linear indexes or the elements above diagonal.
%
%             If input argument A is a scalar, y is a vector of length
%             A*(A-1)/2 containing the (A*(A-1)/2) linear indexes of the
%             elements above diagonal of a square matrix of order A.
%
%             If input argument A is a square matrix of order k (k>1), y
%             is a vector of length (k*(k-1)/2) containing the elements
%             above diagonal of matrix A.
%
%       Remark: if the elements abova diagonal must be extracted
%           repeatedly, it is convenient to find the linear indexes of the
%           elements above diagonal once and for all
%
% Copyright 2008-2014.
% Written by FSDA team
%
%
%<a href="matlab: docsearchFS('upmat2vec')">Link to the help function</a>
% Last modified 08-Dec-2013

% Examples:

%{
    % Extract the elements above the main diagonal of a square matrix
    % of order 5
    A=magic(5);
    disp('Input matrix A')
    disp(A)
    disp('Extract the elements above diagonal of matrix A in a vector')
    y=Upmat2vec(A);
    disp(y)
%}

%{
    % Finds the linear indexes of the elements above the main diagonal of a
    % square matrix of order 5
    k=5;
    disp(['Input scalar is k=' num2str(k)])
    disp('Extract the linear indexes of the elements above diagonal')
    disp(['of a square matrix of size ' num2str(k) ' in a vector'])
    y=Upmat2vec(k);
    disp(['The linear indexes of the elements above diagonal of a square matrix of size ' num2str(k) ' are:'])
    disp(y)
    disp(['For example if A=invhilb(' num2str(k) ')'])
    A=invhilb(k)
    disp('The elements above diagonal of A are');
    A1=A(:);
    disp(A1(y))
%}


%% Beginning of code

if isscalar(A)
    k=A;
elseif ismatrix(A)
    [k,p]=size(A);
    if k ~= p
        error('Input matrix must be square')
    end
else
    error('Input matrix must be a scalar or a square matrix')
end

ind=(1:(k*(k-1)/2))';

j = round(floor(-.5 + .5 * sqrt(1 + 8 * (ind - 1))) + 2);
i = round(j .* (3 - j) / 2 + ind - 1);
% y contains the linear indexes of the upper tringualar part of a square
% matrix of size k
y = i + (j - 1).*k;


if ~isscalar(A)
    % Extract the elements above diagonal of matrix A
    A1=A(:);
    y=A1(y);
end

end
