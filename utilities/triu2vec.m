function y=triu2vec(A,k)
%triu2vec extracts in a vector the linear indexes or the elements on and above the k-th diagonal of a square matrix
%
% k is an optional non-negtive integer that, if not given, is set to 0 to
% indicate the main diagonal. The indices are extracted following the
% traditional "packed storage" scheme for symmetric, Hermitian or
% triangular matrices, adopted by MATLAB for the linear indexing.
%
%<a href="matlab: docsearchFS('triu2vec')">Link to the help function</a>
%
%  Required input arguments:
%
%         A   : scalar, or square matrix  of order r (r>1). Scalar or
%               square matrix.
%               If A is a non negative integer, it indicates the number of
%               rows (columns) of a matrix whose linear indexes of the
%               elements on and above diagonal have to be found. Then,
%               triu2vec(A) returns the (A*(A-1)/2)+A linear indexes of the
%               elements on and above diagonal.
%
%               If A is a square matrix of order r (r>1), then triu2vec(A)
%               returns the elements on and above diagonal of matrix A in a
%               vector of length (r*(r-1)/2)+r;
%
%  Optional input arguments:
%
%         k   : which diagonal. Scalar. Non negative integer given to return the elements on and
%               above the k-th diagonal of A, being k=0 the main diagonal.
%               Default is k = 0, i.e., the main diagonal is also returned.
%               Negative integers are treated as 0, i.e. elements on and
%               above the main diagonal are returned. No linear indices or
%               elements are returned if the user provides an integer k
%               larger than the order of the input matrix.
%                 Example - 'k',0
%                 Data Types - double
%
%  Output:
%
%        y  : vector containing the linear indexes or the elements on and
%             above the k-th diagonal.
%
%                     For example if k =0 or nargin ==1
%             If input argument A is a scalar, y is a vector of length
%             (A*(A-1)/2)+A containing the linear indexes of the
%             elements on and above diagonal of a square matrix of order A.
%
%             If input argument A is a square matrix of order r (r>1), y
%             is a vector of length (r*(r-1)/2)+r containing the elements
%             on and above diagonal of matrix A.
%
%                   For example if k == 2
%             If input argument A is a scalar, y is a vector of length
%             (A*(A-1)/2)+A -A -(A-1) containing the linear indexes of the
%             elements on and above the 2nd diagonal of a square matrix of
%             order A.
%
%             If input argument A is a square matrix of order r (r>1), y
%             is a vector of length (r*(r-1)/2)+r -r -(r-1) containing the
%             elements on and above the 2nd diagonal of a square matrix of
%             order A.
%
%
%       Remark 1: if the elements on and above diagonal must be extracted
%             repeatedly, it is convenient to find the linear indexes once
%             and for all
%
%       Remark 2: the linear indices of the upper triangualar part of a
%       square matrix A of size r can be also obtained with the following
%       formulas and code:
%
%               ind=(1:(r*(r-1)/2))';
% 
%                 j = round(floor(-.5 + .5 * sqrt(1 + 8 * (ind - 1))) + 2);
%                 i = round(j .* (3 - j) / 2 + ind - 1);
% 
%                 y = i + (j - 1).*r; % the linear indexes
% 
%                 A1 = A(:);
%                 y  = A1(y);         % the elements above diagonal of matrix A
%
%
%
% See also: diag
%
% Copyright 2008-2016.
% Written by FSDA team
%
%<a href="matlab: docsearchFS('triu2vec')">Link to the help function</a>
% Last modified lun 16 mag 2016 23:43:20

% Examples:

%{
    % Extract the elements above the main diagonal. 
    % The input is a square matrix of order 5.
    A=magic(5);
    disp('Input matrix A')
    disp(A)
    disp('Extract the elements above the main diagonal of matrix A in a vector')
    k = 1;
    y=triu2vec(A,k);
    disp(y)
%}

%{
    % Extract the elements on and above the main diagonal. 
    % The input is a square matrix of order 5.
    A=magic(5);
    disp('Input matrix A')
    disp(A)
    disp('Extract the elements on and above diagonal of matrix A in a vector')
    y=triu2vec(A);
    disp(y)
%}

%{
    % Finds the linear indexes of the elements on and above the main diagonal of a
    % square matrix of order 5
    r=5;
    disp(['Input scalar is r=' num2str(r)])
    disp('Extract the linear indexes of the elements on and above diagonal')
    disp(['of a square matrix of size ' num2str(r) ' in a vector'])
    y=triu2vec(r);
    disp(['The linear indexes of the elements on and above diagonal of a square matrix of size ' num2str(r) ' are:'])
    disp(y)
    disp(['For example if A=invhilb(' num2str(r) ')'])
    A=invhilb(r)
    disp('The elements on and above diagonal of A are');
    A1=A(:);
    disp(A1(y))
%}


%% Beginning of code

if isscalar(A) && rem(A,1) == 0
    r=A;
    sizeA=[r, r];
elseif ismatrix(A) && ~isscalar(A) 
    [r,c]=size(A);
    sizeA=[r,c];
    if r ~= c
        error('FSDA:triu2vec:WrongInput','Input matrix must be square')
    end
else
    error('FSDA:triu2vec:WrongInput','Input must be a non negative integer or a square matrix')
end

if nargin == 1 || rem(k,1) ~= 0 || k<0
    k=0;
end

maskA	= logical(triu(ones(sizeA),k));
if isscalar(A)
    y = find(maskA);
else
    y	= A(maskA);
end

if isempty(y)
    disp(['Nothing to return on and above the ' num2str(k) 'th diagonal, for a matrix of order ' num2str(r) ' (being 0 the main diagonal)']);
end

end
%FScategory:UTIGEN