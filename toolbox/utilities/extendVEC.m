function [ext] = extendVEC(x, index, value)
% extendVEC extends vector x inserting value in the positions in index
%
%<a href="matlab: docsearchFS('extendVEC')">Link to the help function</a>
%
% Required input arguments:
%
% x:        Vector to be extended. Numeric vector.
%           It can be either a row or a column vector.
%
% index:    A vector containing the positions where to add value inside x.
%           Numeric vector with integer elements.
%           The index values are sorted so that the code starts to insert
%           the element "value" from the smallest position.
%
% value:    The element to insert in x for the extension. Scalar. 
%           Any numeric scalar, including NaN, Inf, and -Inf is allowed
%
%
% Optional input arguments:
%
% Output:
%
% ext:  The extension of vector x with the element value inserted in the positions in index.
%       For index > length(x), the code inserts index-length(x) elements
%       equal to value in the tail of x. It can be either a row or a column
%       vector, depending on input x.
%
%
% Copyright 2008-2025.
% Written by FSDA team
%
%
% See also:
%
%<a href="matlab: docsearchFS('extendVEC')">Link to the help function</a>
%
%$LastChangedDate:: 2020-06-09 16:58:50 #$: Date of the last commit
%
% Example:

%{
  %% Row vector extended with the inclusion of multiple elements "Inf" in the tail.
  x = 1:6;
  index = [2 5 11];
  ext = extendVEC(x, index, Inf);
%}

%{
% Column vector extended with the inclusion of multiple elements "Inf" in the tail.
  x = (1:6);
  index = [2 5 11];
  ext = extendVEC(x, index, Inf);
%}


%{
  % vector extended with the inclusion of multiple elements 99 in the tail.
  % In this case index is a column vector.
  x = (1:6);
  index = [3 5 6]';
  ext = extendVEC(x, index, 99);
%}

%% Beginning of code

% The first argument which is passed is x.
if nargin<1 || isempty(x)==1
    error('FSDA:extendVEC:MissingInputs','Input vector x not specified.');
end
validateattributes(x,{'numeric'},{'vector'})

% The second argument which is passed is index.
if nargin<2 || isempty(index)==1
    error('FSDA:extendVEC:MissingInputs','Input vector index not specified.');
end
% index = numeric vector with integer elements
validateattributes(index, {'numeric'}, {'vector', 'integer'});

% The third argument which is passed is value.
if nargin<3 || isempty(value)==1
    error('FSDA:extendVEC:MissingInputs','Input value not specified.');
end

% Any numeric scalar, including NaN, Inf, and -Inf is allowed
validateattributes(value, {'numeric'}, {'scalar'});


index = sort(index(:))'; % Make sure that index becomes a row vector
ext = x(:)'; % Make sure theat ext becomes a row vector

for r = index
    if r == 1
        ext = [value ext]; %#ok<*AGROW>
    elseif r <= length(ext)
        ext = [ext(1:(r-1)) value ext(r:end)];
    else
        ext = [ext repmat(value, 1, r-length(ext))];
    end
end
% Force ext to be either a row or column vector depending on the initial
% dimension of x
if size(x,1)>1
    ext=ext';
end

end